function t_satrange % ‰q¯-ƒ‚ƒjƒ^‹ÇŠÔ‹——£’P‘ÌŒ±

C=299792458;
t=caltomjd([2008,1,1,6,0,0]);
state=eletostate([42164000,0.1,45,327,270,90]);
rr=GeodToEcef([40,135,0]);
dtr=0.0001;
U=EcsfToEcef(t);

% LightTime•â³OFF
disp('satrange:corrlightt=0')
[rs,range,drds]=satrange(state,rr,U,dtr,0);
rs,range,drds
%range_r=norm(U*state(1:3)-rr);
%drds_r=drds_diff(state,rr,U,dtr,0);
%range-range_r, range./range_r
%drds-drds_r, drds./drds_r


% LightTime•â³ON
disp('satrange:corrlightt=1')
[rs,range,drds]=satrange(state,rr,U,dtr,1);
rs,range,drds
%drds_r=drds_diff(state,rr,U,dtr,1);

%range_r=norm(rs-EcsfToEcef(t-dtr/86400)'*rr);
%drds_r=drds_diff(state,rr,U,dtr,1);
%drdt_r=drdt_diff(state,rr,U,dtr,1);
%range-range_r, range./range_r
%drds-drds_r, drds./drds_r
%drdt_r

% ·•ª‹ß— --------------------------------------------------------------------
function drds=drds_diff(state,rr,U,dtr,corrlight)
[rs,range]=satrange(state,rr,U,dtr,corrlight);
ds=diag([0.01,0.01,0.01,0.0001,0.0001,0.0001]);
for n=1:6
    [rs,rng]=satrange(state+ds(:,n),rr,U,dtr,corrlight);
    drds(n)=(rng-range)/ds(n,n);
end

function drdt=drdt_diff(state,rr,U,dtr,corrlight)
C=299792458;
[rs,range]=satrange(state,rr,U,dtr,corrlight);
[rs,rng]=satrange(state,rr,U,dtr+0.00001,corrlight);
drdt=(rng-range)/0.00001/C; % dƒÏ/dcdtr=dƒÏ/ddtr/c
