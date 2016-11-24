function gpsrtk(prm)
%-------------------------------------------------------------------------------
% [system] : GpsTools
% [module] : RTK-GPS estimator
% [func]   : estimate station position by RTK(real-time kinematic)
% [func]   : detect/repair cycle-slip, smooth codes and generate clean obs. data
% [argin]  : (prm) = processing parameters struct (see source-code)
% [argout] : none
% [note]   : output files
%            posf_{rcv}_YYYYMMDDHH.mat : receiver position
%              epoch = start epoch [year,month,day,hour,min,sec])
%              time  = relative time vector to epoch (sec)
%              rcv   = receiver name
%              data  = estimated states (data(n,:) : time(n) estimation)
%              covs  = varience of est. (covs(n,:) : time(n) varience)
%              prm   = processing parameters struct
%
% [version]: $Revision: 12 $ $Date: 2008-11-25 10:02:15 +0900 (火, 25 11 2008) $
%            Copyright(c) 2004-2006 by T.Takasu, all rights reserved
% [history]: 05/12/24  0.1  new
%-------------------------------------------------------------------------------
if nargin<1, prm=[]; end, p=setprm(prm);

% read observation data
[zz,iz]=readobs(p.t0,p.t,p.sats,p.rcvs,p.dir.obs);
if isempty(zz), error('observation data empty'), end

% read navigation messages
[nav,inav]=readnav(p.t0,p.t,p.sats,{},p.dir.nav);
if isempty(nav), error('navigation messages empty'), end

% read reference station positions
rref=readpos(p.t0,p.t(1),p.rcvs(2:end),p.dir.pos,p.pos);

% estimate states
[xs,ps,xf,pf,xp,ix]=eststate(zz,iz,nav,inav,rref,p);

% save states
savestate(xf,ps,ix,p);

% plot results
if p.plot(1), plot_pos('(FLOAT)',xs,ps,rref,p,1); end
if p.plot(2), plot_pos('(FIX)',xf,pf,rref,p,1); end
if p.plot(3), plot_pos('(PP)',xp,ps,rref,p,0); end
if p.plot(4), plot_ion(xf,ps,ix,p); end

% set default parameters -------------------------------------------------------
function p=setprm(prm)

% default parameters
p.t0=caltomjd([2004,10,1]); % day (mjd-gpst)
p.t    =0:30:3570;          % time of day (sec)
p.sats ={};                 % satellites
for n=1:31, p.sats={p.sats{:},sprintf('GPS%02d',n)}; end
p.rcvs ={'960759','93040'}; % stations (rcvs{1}:positioning sta)
p.mode =1;                  % positioning mode (0=static,1=kinematic,2=dynamic)
p.obs  =1:4;                % observables (1=L1,2=L2,1:2=L1+L2,1:4=L1+L2+P1+P2)
p.trop =0;                  % tropos delay (0=off,1=estimate)
p.iono =0;                  % ionos delay (0=off,1=estimate)
p.pos  ='gsipos';           % reference position ('gsipos'=GSI,'igssnx'=IGS)
p.eph  ='';                 % satellite ephemeris ('':broadcast,'igp':igs-pred)
p.plot =0;                  % plot results (0=off,1=on)
p.ambf =2;                  % ambiguity fixing(0=float,1=round,2=lambda)
p.ambv =2;                  % ambiguity validation(0=no,1=ratio,2=f-ratio,3=chi2)
p.ambc =2;                  % ambiguity validation critical value
p.ambp =1;                  % ambiguity partical fix(0=inst,1=all,2=partal)
p.outp =3;                  % outlier threshold of prefit-residuals (sigmas)
p.outf =3;                  % outlier threshold of postfit-residuals (sigmas)
p.elw  =0;                  % measurement noise elevation weighting(0=off,1=on)
p.elmin=15;                 % elevation cutoff of valid observation (deg)
p.elamb=15;                 % elevation cutoff of ambiguity resolution (deg)
p.slipc=5;                  % cycle-slip counter threshold
p.sig=[0.003,0.3];
p.p=[10,1,0.3,0.1,20];      % initial state/proess noise std. devs
p.q=[0,0,1E-4,1E-3,0];      %    [pos(m),vel(m/s),trop(m),ion(m),amb(cycle)]
p.f=[1.57542E9,1.2276E9];   % carrier frequencies (Hz)
p.plot=[0,0,0,0];           % plot result flags
p.range=[-0.1,0.1];         % plot y-axis range

p.dir.obs='l:\gps\obs_gsi\%Y%m%d'; % observation data directory
p.dir.nav='k:\gps\nav\%Y%m';       % navigation messages directory
p.dir.pos='k:\gps\pos';            % position data directory
p.dir.eph='k:\gps\sp3\%Y%m';       % satellite ephemeris directory
p.dir.est='';                      % output data directory

for f=fieldnames(p)'
    if isfield(prm,f{1}), p=setfield(p,f{1},getfield(prm,f{1})); end
end

b2R={'GPS02','GPS11','GPS13','GPS14','GPS16','GPS18','GPS19','GPS20','GPS21',...
     'GPS22','GPS23','GPS28'};
for n=1:length(p.sats), p.satt(n)=any(strcmp(p.sats{n},b2R)); end

% estimate states --------------------------------------------------------------
function [xs,ps,xf,pf,xp,ix]=eststate(obs,iobs,nav,inav,rref,p)

[time,ti]=intersect(round(iobs(:,1)),p.t); ti=[0,ti];

[x,P,ix]=initstate(p);

ref=[]; slip=zeros(length(ix.n),1);

xs=repmat(nan,length(time),length(x)); ps=xs; xf=xs; pf=xs;
xp=zeros(length(time),3);

for k=1:length(time)
    i=ti(k)+1:ti(k+1); obsk=obs(i,:); iobsk=iobs(i,:);
    
    % approx receiver position/clock-bias/elevation
    [rr,dt,el]=pointp(time(k),obsk,iobsk,nav,inav,p);
    
    % select reference satellie
    [x,P,ref]=selref(x,P,ix,ref,obsk,iobsk,el,p);
    
    % temporal update of states
    [x,P,is,va]=udstate(time(k),x,P,ix,rr,obsk,iobsk,el,ref,slip,p);
    
    % double-difference of observables
    z=obs_dd(time(k),obsk,iobsk,p.obs,el,ref,p);
    
    % measurement model by prefit states
    [h,H,R]=measmodel(time(k),x,ix,rref,dt,nav,inav,el,ref,p);
    
    for outl=0:length(z)-3
        % measurement update of states
        [xu,Pu]=filt(x,P,z,h,H,R,p);
        
        % postfit residuals
        s=z-measmodel(time(k),xu,ix,rref,dt,nav,inav,el,ref,p);
        
        % exclude outliers
        [maxs,i]=max(s.^2./diag(R));
        if p.outf==0|maxs<p.outf^2, x=xu; P=Pu; break, else z(i)=nan; end
        
        % increment slip-counter of invalid phase-obs
        if is(i)>0, slip(is(i))=slip(is(i))+1; end
    end
    % decrement slip-counter of valid phase-obs
    for i=find(is(:)>0&~isnan(z))', slip(is(i))=max(slip(is(i))-1,0); end
    
    xs(k,:)=x'; ps(k,:)=diag(P)'; xp(k,:)=rr(:,1)';
    
    % fix ambiguity
    if ~isempty(h)
        i=find(~isnan(s)); res=s(i)'/R(i,i)*s(i);
        j=ix.n; j=j(~isnan(x(j))&x(j)~=0&repmat(va,ix.nm,1)&slip<=0);
        [xx,Px,f,fix]=fixamb(x,ix,P,res,j,p);
        if fix, xf(k,:)=xx'; pf(k,:)=diag(Px)'; end
        if fix&p.ambp==2
            i=find(~isnan(xx)); x(i)=xx(i); P(:,i)=Px(:,i); P(i,:)=Px(i,:);
        end
    end
    e=mjdtocal(p.t0,time(k)); fs='* ';
    msg('%02d:%02d:%02.0f: %12.3f %12.3f %12.3f%c: r=%s n=%3d,%3d s=%5.1f f=%5.1f %5.1f',...
        e(4:6),xx(ix.r),fs(fix+1),p.sats{ref},length(z(~isnan(z))),outl,res,f);
end

% initialize states ------------------------------------------------------------
function [x,P,ix]=initstate(p)
n=length(p.rcvs); m=length(p.sats)*(n-1); k=length(find(intersect(p.obs,1:2)));
ix.r=1:3; nx=3;
if p.mode==2, ix.v=nx+(1:3); nx=nx+3; else ix.v=[]; end
if p.trop==1, ix.t=nx+(1:n); nx=nx+n; else ix.t=[]; end
if p.iono==1, ix.i=nx+(1:m); nx=nx+m; else ix.i=[]; end
ix.n=nx+(1:k*m); nx=nx+k*m; ix.nn=m; ix.nm=k;
x=repmat(nan,nx,1); P=zeros(nx);

% approx receiver position/clock-bias/elevation by point positiongin -----------
function [rr,dt,el]=pointp(t,obs,iobs,nav,inav,p)
C=299792458; ff=[p.f(1)^2;-p.f(2)^2]/(p.f(1)^2-p.f(2)^2);
for m=1:length(p.sats), rs(:,m)=navtostate(p.t0,t,nav(inav==m,:)); end
for n=1:length(p.rcvs)
    i=find(iobs(:,3)==n&all(~isnan(obs(:,3:4)),2));
    [rr(:,n),dt(n)]=pointpos(p.t0,t,obs(i,3:4)*ff,iobs(i,2:3),nav,inav);
    if ~isempty(i), dt(n)=dt(n)/C-(iobs(i(1),1)-t); end
    for m=1:length(p.sats)
        azel=satazel(rs(:,m),rr(:,n)); el(n,m)=azel(2)*180/pi;
    end
end

% select reference satellite ---------------------------------------------------
function [x,P,ref]=selref(x,P,ix,ref,obs,iobs,el,p)
for n=1:length(p.sats) % obs counts
    nz(n)=length(find(iobs(:,2)==n&all(~isnan(obs(:,p.obs)),2)));
end
% max elevation satellite
i=find(nz>=length(p.rcvs)); [maxel,j]=max(el(1,i)); new=i(j);
if isempty(new)|(~isempty(ref)&ref==new), return, end

A=eye(length(x)); j=0; k=0; nn=length(p.sats);
for n=2:length(p.rcvs)
    if p.iono==1
        i=ix.i(j+(1:nn)); x(i(ref))=0; A(i,i(new))=A(i,i(new))-1; j=j+nn;
    end
    for m=intersect(p.obs,1:2)
        i=ix.n(k+(1:nn)); x(i(ref))=0; A(i,i(new))=A(i,i(new))-1; k=k+nn;
    end
end
i=find(~isnan(x)); x(i)=A(i,i)*x(i); P(i,i)=A(i,i)*P(i,i)*A(i,i)'; ref=new;

% temporal update of states ----------------------------------------------------
function [x,P,is,va]=udstate(t,x,P,ix,rr,obs,iobs,el,ref,slip,p)
C=299792458; vn=[]; va=[]; s=[];

% valid signle-diff/double-diff/ambiguity
for n=2:length(p.rcvs)
    for m=1:length(p.sats)
        vn=[vn;m~=ref&all(el([1,n],m)>=p.elmin)];
        va=[va;m~=ref&all(el([1,n],m)>=p.elamb)];
        if all(el([1,n],m)>=p.elmin), s=[s;(n-2)*length(p.sats)+m]; end
    end
end
is=[]; for i=p.obs, if i<=2, is=[is;s]; s=s+ix.nn; else is=[is;s*0]; end, end

F=eye(length(x));
Q(ix.r)=p.q(1); Q(ix.v)=p.q(2); Q(ix.t)=p.q(3); Q(ix.i)=p.q(4); Q(ix.n)=p.q(5);

if p.mode==1|any(isnan(x(ix.r)))
    x(ix.r)=rr(:,1); P(:,ix.r)=0; P(ix.r,:)=0; Q(ix.r)=p.p(1);
end
if p.mode==2
    if any(isnan(x(ix.v))), x(ix.v)=0; P(:,ix.v)=0; P(ix.v,:)=0; Q(ix.v)=p.p(2);
    else x(ix.r)=x(ix.r)+x(ix.v)*(p.t(2)-p.t(1)); end
end
if p.trop==1
    i=ix.t(isnan(x(ix.t))); x(i)=0; P(:,i)=0; P(i,:)=0; Q(i)=p.p(3);
end
if p.iono==1
    i=ix.i(vn&isnan(x(ix.i))); x(i)=0; P(:,i)=0; P(i,:)=0; Q(i)=p.p(4);
    i=ix.i(~vn); x(i)=nan; F(i,i)=0; Q(i)=0;
end
if p.ambp==0, x(ix.n)=nan; P(:,ix.n)=0; P(ix.n,:)=0; Q(ix.n)=0; end
ch=intersect(1:2,p.obs); fs=[];
for n=2:length(p.rcvs)
    for m=intersect(1:2,p.obs), fs=[fs;ones(length(p.sats),1)*p.f(m)/C]; end
end
i=find(repmat(vn,ix.nm,1)&isnan(x(ix.n)));
if ~isempty(i)
    N=(obs_dd(t,obs,iobs,ch,[],ref,p)-obs_dd(t,obs,iobs,ch+2,[],ref,p)).*fs;
    x(ix.n(i))=N(i); P(:,ix.n(i))=0; P(ix.n(i),:)=0; Q(ix.n(i))=p.p(5);
    msg('reinit ambiguity :%s',sprintf(' %d',i));
end
i=ix.n(~repmat(vn,ix.nm,1)); x(i)=nan; P(:,i)=0; P(i,:)=0; Q(i)=0;

P=P+diag(Q.^2);

% measurement update of states -------------------------------------------------
function [x,P]=filt(x,P,y,h,H,R,p)
if isempty(y), return, end
s=y-h;
i=find(~isnan(s)); j=find(~isnan(x)&diag(P)>0);
s=s(i); H=H(i,j); S=H*P(j,j)*H'+R(i,i);
if p.outp>0 % prefit residuals test
    k=find(s.^2<diag(S)*p.outp);
    s=s(k); H=H(k,:); S=S(k,k);
end
K=P(j,j)*H'/S;
x(j)=x(j)+K*s;
P(j,j)=P(j,j)-K*H*P(j,j);

% double difference of observables ---------------------------------------------
function y=obs_dd(t,obs,iobs,ch,el,ref,p)
y=[];
for n=2:length(p.rcvs)
    yn=[];
    y11=obsdat(t,ref,1,obs,iobs,ch,p);
    y21=obsdat(t,ref,n,obs,iobs,ch,p);
    for m=1:length(p.sats)
        if isempty(el)|all(el([1,n],m)>=p.elmin)
            y12=obsdat(t,m,1,obs,iobs,ch,p);
            y22=obsdat(t,m,n,obs,iobs,ch,p);
            yn=[yn;(y11-y12)-(y21-y22)];
        end
    end
    y=[y;reshape(yn,prod(size(yn)),1)];
end

% retrieve observation data ----------------------------------------------------
function z=obsdat(t,sat,rcv,zz,iz,obs,p)
C=299792458;
i=find(iz(:,2)==sat&iz(:,3)==rcv);
if ~isempty(i), z=zz(i,1:4).*[C./p.f,1,1]; else z=[nan,nan,nan,nan]; end
z=z(obs);

% fix ambiguity ----------------------------------------------------------------
function [x,P,f,fix]=fixamb(x,ix,P,s,j,p);
i=[ix.r,ix.i]; i=i(~isnan(x(i))); r=0; f=[0,0]; fix=0;
switch p.ambf % ambiguity fixing strategy
case 0, return                              % float
case 1, N=round(x(j)); f=[1,1000];          % round
%case 2, [N,f]=lambda(x(j),P(j,j),2);        % lambda
case 2, [N,f]=lambda(x(j),P(j,j),2,1);        % lambda
end
n=length(j); c=chidist(n,0.005);
switch p.ambv % ambiguity validation
case 0, fix=1;                              % no validation
case 1, fix=f(2)/f(1)>=p.ambc;              % ratio test
case 2, fix=(s+f(2))/(s+f(1))>=p.ambc;      % f-ratio test
case 3, fix=f(1)/n<c;
case 4, fix=f(1)/n<c&f(2)/n>c;
end
if fix
    x([i,j])=[x(i)-P(i,j)/P(j,j)*(x(j)-N(:,1));N(:,1)];
    P(i,i)=P(i,i)-P(i,j)/P(j,j)*P(j,i); P(:,j)=0; P(j,:)=0;
end

% phase/code measurements model ------------------------------------------------
function [h,H,R]=measmodel(t,x,ix,rref,dt,nav,inav,el,ref,p)
C=299792458; lam=C./p.f; fi=[1,p.f(1)^2/p.f(2)^2]; fi=[fi,-fi];
h=[]; H=[]; R=[]; i=0; j=0; nn=length(p.sats);
for n=2:length(p.rcvs)
    [r,dr,Rr,k]=...
        geo_dd(t,[x(ix.r),rref(:,:,n-1)'],dt([1,n]),nav,inav,el([1,n],:),ref,p);
    for m=p.obs
        Hi=[dr,zeros(size(dr,1),length([ix.t,ix.i,ix.n]))];
        if m<=2 % phase
            hi=r+lam(m)*x(ix.n(i+k)); Ri=Rr*p.sig(1)^2;
            Hi(:,ix.n(i+k))=lam(m)*eye(length(k)); i=i+nn;
        else % code
            hi=r; Ri=Rr*p.sig(2)^2;
        end 
        if p.iono==1
            hi=hi+fi(m)*x(ix.i(j+k));
            Hi(:,ix.i(j+k))=fi(m)*eye(length(k));
        end
        if ~isempty(hi), h=[h;hi]; H=[H;Hi]; R=blkdiag(R,Ri); end
    end
    j=j+nn;
end

% double-difference of geometry ------------------------------------------------
function [r,dr,R,i]=geo_dd(t,rr,dt,nav,inav,el,ref,p)
nav1=nav(inav==ref,:);
[r11,dr11]=geoterm(t-dt(1),rr(:,1),ref,nav1,el(1,ref),p);
[r21,dr21]=geoterm(t-dt(2),rr(:,2),ref,nav1,el(2,ref),p);
w1=1/sin(el(1,ref))^2+1/sin(el(2,ref))^2;
r=[]; dr=[]; w=[]; i=[];
for m=1:length(p.sats)
    if all(el(1:2,m)>=p.elmin)
        nav2=nav(inav==m,:);
        [r12,dr12]=geoterm(t-dt(1),rr(:,1),m,nav2,el(1,m),p);
        [r22,dr22]=geoterm(t-dt(2),rr(:,2),m,nav2,el(2,m),p);
        r=[r;dd(r11,r12,r21,r22)]; dr=[dr;dr11-dr12];
        w=[w;1/sin(el(1,m))^2+1/sin(el(2,m))^2]; i=[i,m];
    end
end
if p.elw, R=ones(length(w))*w1+diag(w);
else R=2*(ones(length(w))+eye(length(w))); end

% geometric term ---------------------------------------------------------------
function [r,dr]=geoterm(t,rr,sat,nav,el,p)
if isempty(p.eph)
    [r,dr]=geodist(p.t0,t,nav,rr);     % by broardcast ephemeris
else
    [r,dr]=geodistp(p.t0,t,sat,rr,p);  % by precise ephemeris
end
r=r+tropos(t,rr,el,p);

% geometric distance by precise-ephemeris --------------------------------------
function [r,dr]=geodistp(t0,t,sat,rr,p)
C=299792458; OMGE=0.7292115E-4; tau=0;
while 1 % solve light-time eq.
    rs=Rz(OMGE*tau)*sateph(t0,t-tau,sat,p);
    r=norm(rr-rs);
    if isnan(r)|abs(r-tau*C)<1E-4, break, else tau=r/C; end
end
dr=(rr-rs)'/r; % partial derivative

% read precise ephemeris -------------------------------------------------------
function rs=sateph(t0,t,sat,p)
persistent te0 te eph i j ti ephi
if isempty(te0)|te0~=p.t0
    te0=p.t0; te=900*(floor(p.t(1)/900)-5:floor(p.t(end)/900)+5);
    eph=readeph(te0,te,p.sats,p.dir.eph,p.eph); ie=[];
end
k=min(find(t<te));
if isempty(i)|k~=i
    i=k; j=i-4:i+4; ti=repmat(te(j)',1,9); tt=repmat(te(j),9,1)-ti;
    for n=1:9, tt(n,n)=1; ephi(n,:,:)=eph(j(n),1:3,:)/prod(tt(:,n)); end
end
tt=t-ti; for n=1:9, tt(n,n)=1; end
rs=(prod(tt,1)*ephi(:,:,sat))';
if ~p.satt(sat), rs=rs-1.023*rs/norm(rs); end % satellite antenna offset

% tropospheric delay -----------------------------------------------------------
function e=tropos(t,rr,el,p);
gpos=eceftogeod(rr); t=p.t0+t/86400;
[zh,zw]=trop_saast(t,[0,pi/2],gpos);     % zenith delays
[mh,mw]=mapf_nmf(t,[0,el*pi/180],gpos);  % mapping functions
e=mh*zh+mw*zw;

% f-distribution ---------------------------------------------------------------
function f=fdist(n,m,a)
r=[0,10]; while a<fp(n,m,r(2)), r(2)=r(2)*4; end
while 1
    f=(r(1)+r(2))/2; p=fp(n,m,f);
    if abs(p-a)<a*1E-3, break, elseif p<a, r(2)=f; else r(1)=f; end
end

function p=fp(n,m,x), p=1-betainc(n*x/(m+n*x),n/2,m/2);

% chi-square distribution ------------------------------------------------------
function x=chidist(n,a)
r=[0,10]; while a<x2p(n,r(2)), r(2)=r(2)*4; end
while 1
    x=(r(1)+r(2))/2; p=x2p(n,x);
    if abs(p-a)<a*1E-3, break, elseif p<a, r(2)=x; else r(1)=x; end
end

function p=x2p(n,x), p=1-gammainc(x/2,n/2);

% save states ------------------------------------------------------------------
function savestate(xs,ps,ix,p)
saves('posf',xs(:,ix.r),ps(:,ix.r),p);
saves('ionf',xs(:,ix.i),ps(:,ix.i),p);

function saves(type,data,covs,p)
epoch=mjdtocal(p.t0,p.t(1)); time=p.t-p.t(1); rcv=p.rcvs{1}; prm=p;
f=fullfile(p.dir.est,sprintf('%s_%s_%04d%02d%02d%02d.mat',type,rcv,epoch(1:4)));
save(f,'epoch','time','rcv','data','covs','prm');

% plot results -----------------------------------------------------------------
function plot_pos(ti,xs,ps,rref,p,flg)
rpos=readpos(p.t0,p.t(1),p.rcvs{1},p.dir.pos,p.pos);
[g,E]=eceftogeod(rpos');
for k=1:size(xs,1)
    errs(k,:)=(xs(k,1:3)-rpos)*E';
    sigs(k,:)=sqrt(diag(E*diag(ps(k,1:3))*E'))';
end
s='';
if flg
    s=[s,'ref='];
    for n=1:length(p.rcvs)-1
        s=[s,sprintf('%s(%.1fkm) ',p.rcvs{n+1},norm(rref(:,:,n)-rpos)/1E3)];
    end
    n=size(xs,1); m=length(find(~isnan(xs(:,1))));
    k=length(find(sum(errs(:,1:2).^2,2)>0.1));
    s=[s,sprintf('fix=%.1f%% mis=%.2f%%',m/n*100,k/n*100)];
end
gut('newfigg','','',[600,400]); yl={'East (m)','North (m)','Up (m)'};
for n=1:3
    if n<3, topts='nolabel'; else topts=''; end
    ggt('subplotv','',n,3,'taxis',p.t0,'xlim',[p.t(1),p.t(end)]/3600,...
        'ylim',p.range,'topts',topts,'margin',[0.08,0.03,0.04,0.012]); 
    plot(p.t/3600,errs(:,n),'.-'); ylabel(yl{n});
%    plot(p.t/3600,-sigs(:,n),'r:'); plot(p.t/3600,sigs(:,n),'r:');
    ggt('mtext',['_',num2str(n)],sprintf('MEAN:%7.4fm RMS:%7.4fm',...
        meann(errs(:,n)),rms(errs(:,n))));
    if n==1
        title(['RTK Position ',ti,' ',p.rcvs{1},' : ',tstr(p.t0),' : ',s]);
    end
end

function plot_ion(xs,ps,ix,p)
gut('newfigg','','',[600,400]);
ggt('subplotv','',1,1,'taxis',p.t0,'xlim',[p.t(1),p.t(end)]/3600,...
    'ylim',[-0.5,0.5],'margin',[0.08,0.03,0.04,0.012]); 
h=plot(p.t/3600,xs(:,ix.i),'.-');
ylabel('fixed solutions (m)'); ylim([-0.1,0.1]); %legend(h,p.sats);
title(['Ionospheric Delay : ',p.rcvs{1},' : ',tstr(p.t0)]);

% double difference  -----------------------------------------------------------
function d=dd(r11,r12,r21,r22), d=(r11-r12)-(r21-r22);

% rotation matrix around z axis ------------------------------------------------
function R=Rz(t), R=[cos(t),sin(t),0;-sin(t),cos(t),0;0,0,1];

% rms/mean ---------------------------------------------------------------------
function r=rms(x), r=sqrt(mean(x(~isnan(x(:,1))).^2,1));
function r=meann(x), r=mean(x(~isnan(x(:,1)),:),1);

% time string ------------------------------------------------------------------
function s=tstr(t), ep=mjdtocal(t); s=sprintf('%04d/%02d/%02d',ep(1:3));

% debug print ------------------------------------------------------------------
function msg(varargin), disp(sprintf(varargin{:}));
