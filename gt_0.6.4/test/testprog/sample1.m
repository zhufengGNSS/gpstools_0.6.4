% çLïÒóÔ/ê∏ñßóÔî‰är
function sample1

dirs.erp='d:\gps\erp';
dirs.eph='d:\gps\eph';
dirs.nav='d:\gps\nav';
dirs.clk='d:\gps\clk';

td=CalToMjd([2003,12,4]);
time=0:900:86400*3;

sats={...
%'GPS01','GPS02','GPS03','GPS04','GPS05','GPS06','GPS07','GPS08','GPS09','GPS10',...
%'GPS11','GPS13','GPS14','GPS15','GPS16','GPS17','GPS18','GPS20','GPS21','GPS23',...
%'GPS24','GPS25','GPS26','GPS27','GPS28','GPS29','GPS30','GPS31'
'GPS29'
};

U=zeros(3,3,length(time));
for n=1:length(time)
    tutc=(td+time(n)+19-32)/86400;
    erp=ReadErp(tutc,dirs.erp)+[ErpVar(tutc,-32),0,0];
    U(:,:,n)=EcsfToEcef(tutc,erp,-32);
end
for n=1:length(sats), PlotEphError(td,time,sats{n},U,dirs), end

function PlotEphError(td,time,sat,U,dirs)
C=299792458;
ephp=ReadEph(td,time,sat,U,dirs.eph);
clkp=C*ReadClk(td,time,sat,dirs.clk);
nav=ReadNav(td,time(1):86400:time(end),sat,{},dirs.nav);
for n=1:length(time)
    [pos,dts]=NavToState(td,time(n),nav);
    ephn(n,:)=(U(:,:,n)'*pos)';
    clkn(n,1)=C*dts(1);
end
deph=repmat(nan,length(time),5);
for n=1:length(time)
    E=EcsfToSatf(ephp(n,:)');
    deph(n,2:4)=(ephn(n,1:3)-ephp(n,1:3))*E';
    deph(n,1)=norm(deph(n,2:4));
    deph(n,5)=clkn(n,1)-clkp(n,1);
end
figure
plot(time/3600,deph)
title(['GPS BROADCAST - IGS PRECISE EPHEMERIS :',sat])
xlabel('time(H)'), ylabel('error(m)')
legend({'total(m)','radial(m)','along-track(m)','cross-track(m)','clock(m)'})
hold on, grid on
disp(sprintf('%s : rms err pos = %.3fm, clk = %.3fm',sat,...
     sqrt(mean(deph(:,1).^2)),sqrt(mean(deph(:,5).^2))))
