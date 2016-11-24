function [x,P,arc,phs,vz,res,iv,ds]=gt_udmeas(td,t,x,ix,P,arc,phs,z,tz,iz,antp,...
                                              state,sig,pass,fb,prm)
%-------------------------------------------------------------------------------
% [system] : GpsTools
% [module] : measurement update of states
% [func]   : measuremnet update of states
% [argin]  : td,t    = date(mjd-gpst), time(sec)
%            x,ix,P  = state variables/state variables index/covariances
%            arc     = arc flags
%            phs     = phase-windup corrections
%            z,tz,iz = obsrevation data/observation data sampling time/index
%            antp    = receiver antenna pcv parameters
%            state   = predetermined states (see gt_readstate.m)
%            pass    = pass number
%            fb      = forward/backward flag(1:forward,2:backward)
%            prm     = processing parameters struct
% [argout] : x,P     = updated state variables/covariances
%            arc     = updated arc flags
%            phs     = updated phase-windup corrections
%            vz      = observation data status (0:invalid,1:valid,2:outlier)
%            res     = residuals
%                      res(n,1)   = iz(n) prefit residual (m)
%                      res(n,2)   = iz(n) postfit residual (m)
%                      res(n,3)   = iz(n) phase bias residual (m)
%            iv      = invalid states index
%            ds      = debug information
%                      ds(n,1:2)  = obs data/model
%                      ds(n,3:5)  = station position
%                      ds(n,6)    = sampling time offset
%                      ds(n,7)    = geometric distance
%                      ds(n,8:9)  = receiver/satellite clock bias
%                      ds(n,10:11)= tropospheric delay dry/wet
%                      ds(n,12)   = phase bias
%                      ds(n,13:14)= satellite/receiver antenna pcv
%                      ds(n,15)   = reletivity correction
%                      ds(n,16)   = phase-windup correction
%                      ds(n,17)   = phase-multipath correction
% [note]   :
% [version]: $Revision: 20 $ $Date: 2009-05-01 04:15:33 +0900 (金, 01 5 2009) $
%            Copyright(c) 2004-2008 by T.Takasu, all rights reserved
% [history]: 06/03/26  0.1  separated from gpsestd.m, marged measmodel.m
%            08/11/24  0.2  support gpt for meterological parameters (gt_0.6.4)
%                           support vmf1 for tropospheric mapping function
%                           fix bug on taking over phase-windup correction
%                           fix bug on receiver attitude model=leo sat
%-------------------------------------------------------------------------------
global utc_tai, utc_tai=prm_utc_tai(td+t/86400,1); tu=td+(t+19+utc_tai)/86400;

% set earth rotation parameter/geocenter offset
if prm.est.erp==1
    i=ix.erp;
    if any(isnan(x(i))), x(i)=state.erp(1:3)'; P(i,i)=diag(sig(i).^2); end
    if ~isempty(i), state.erp(1:3)=x(i); end
    if prm.erpvar, state.erp(1:3)=state.erp(1:3)+erpvar(tu,utc_tai); end
end
i=ix.eco; if any(isnan(x(i))), x(i)=0; P(i,i)=diag(sig(i).^2); end

% eci to ecef tranformation matrix
[state.U,q,q,state.gmst,state.dx,state.dy,state.du]=...
    ecsftoecef(tu,state.erp,utc_tai,prm.nutmodel);

% transform satellite ecef position/velocity to eci
if prm.est.erp==1
    [state.ephs,state.ephv]=toeci(state.ephs,state.ephv,state.U,state.du);
end
% set reference clock
i=find(strcmp(prm.rcvs,prm.clkref));
if isempty(i), state.clkref=0; else state.clkref=state.clkr(1,i); end

if strcmp(prm.mapf,'mapf_vmf1'), state.mfcr=repmat(nan,2,length(prm.rcvs)); end

% geodetic positions
for n=1:length(prm.rcvs)
    gpos(n,:)=eceftogeod(state.posr(:,n));
end
% meteorological parameters
if isempty(state.metr)
    for n=1:length(prm.rcvs)
        if strcmp(prm.metsrc,'gpt')
            [pres,temp,undu]=trop_gpt(tu,gpos(n,:));
            state.metr(:,n)=[pres;temp;prm.metprm(3)];
        else
            state.metr(:,n)=[prm.metprm(1)*(1-2.2557E-5*gpos(n,3))^5.2568;...
                             prm.metprm(2)-6.5E-3*gpos(n,3);prm.metprm(3)];
        end
    end
end
% vmf coefficients
if strcmp(prm.mapf,'mapf_vmf1')
    for n=1:length(prm.rcvs)
        [state.mfcr(1,n),state.mfcr(2,n)]=readvmf(tu,gpos(n,:),prm.dirs.trop);
    end
end
% station displacements
psun=state.U*state.rsun; pmoon=state.U*state.rmoon;
for n=1:length(prm.rcvs)
    state.dpos(:,n)=sitedisp(tu,state.posr(1:3,n),psun,pmoon,prm.odisp(:,:,n),...
                             prm.ophas(:,:,n),state.gmst,state.erp,prm.sitedisp);
end
% reset satellite/station states
izv=iz(~isnan(z(:,1)),:);
[x,P]=resetstats(td,t,x,ix,P,izv,state,sig,prm);
[x,P]=resetstatr(td,t,x,ix,P,izv,state,sig,prm);

% reset phase bias of new arc
for i=find(arc(:,fb)==t)'
    j=ix.arcn{arc(i,3),arc(i,4)}; phs(arc(i,3),arc(i,4))=arc(i,9);
    s=arc(i,8); if s==0, s=sig(j); end
    if ~isempty(j), x(j)=arc(i,7); P(j,:)=0; P(:,j)=0; P(j,j)=s^2; end
end

% filter
[x,P,phn,vz,res,ds]=filt(td,t,x,ix,P,z,tz,iz,antp,state,phs,pass,prm);
i=find(~isnan(phn)); phs(i)=phn(i);

% record phase bias of arc end
iv=[];
for i=find(arc(:,3-fb)==t)'
    j=ix.arcn{arc(i,3),arc(i,4)}; bcp=0; s=0; ph=phs(arc(i,3),arc(i,4));
    if ~isempty(j), bcp=x(j); s=sqrt(P(j,j)); iv=[iv,j]; end
    arc(i,7:9)=[bcp,s,ph];
end
% set nan to non-updated clock/position states
j=[ix.satc{:}]; if ~isempty(j)&prm.sclkmodel==1, j=j(1:2:end); end
k=[ix.rcvc{:}]; if ~isempty(k)&prm.rclkmodel==1, k=k(1:2:end); end
i=[j,k,ix.rcvp{:}]; i=i(diag(P(i,i))>=sig(i)'.^2); x(i)=nan; P(i,i)=nan;

% transform satellite ecef position/velocity to eci ----------------------------
function [ephs,ephv]=toeci(ephs,ephv,U,du)
ephs(4:6,:)=U'*ephs(4:6,:)+du'*ephs(1:3,:);
ephs(1:3,:)=U'*ephs(1:3,:);
if ~isempty(ephv)
    ephv(1:3,:)=(U'*sqrt(ephv(1:3,:))).^2;
    ephv(4:6,:)=(U'*sqrt(ephv(4:6,:))).^2;
end

% reset satellite states -------------------------------------------------------
function [x,P]=resetstats(td,t,x,ix,P,iz,state,sig,prm)
if isempty([ix.sata{:}]), return, end
C=299792458; i=[];
for n=1:length(prm.sats)
    j=ix.sato{n}; k=ix.sats{n}; l=ix.satc{n};
    if ~isempty(j)&isnan(x(j(1)))
        x(j)=state.ephs(:,n); i=[i,j];
    end
    if ~isempty(k)&isnan(x(k(1)))
        x(k)=initobtp(n,prm); i=[i,k];
    end
    if ~isempty(l)&isnan(x(l(1)))
        x(l)=0; x(l(1))=(state.clks(1,n)-state.clkref)*C; i=[i,l]; 
    end
end
if ~isempty(i)
    j=find(~isnan(x)); P(i,j)=0; P(j,i)=0; k=sub2ind(size(P),i,i); P(k)=sig(i).^2;
end

% reset receiver states --------------------------------------------------------
function [x,P]=resetstatr(td,t,x,ix,P,iz,state,sig,prm)
if isempty([ix.rcva{:}]), return, end
C=299792458; i=[];
for n=1:length(prm.rcvs)
    j=iz(iz(:,2)==n,1);
    if sum(~isnan(state.clks(1,j)))>=prm.minobs
        j=ix.rcvp{n}; k=ix.rcvc{n}; l=ix.rcvz{n}; m=ix.rcvg{n};
        if ~isempty(j)&isnan(x(j(1)))
            x(j)=initrcvp(n,state,prm); i=[i,j];
        end
        if ~isempty(k)&isnan(x(k(1)))
            x(k)=0; x(k(1))=(state.clkr(1,n)-state.clkref+state.clko)*C; i=[i,k]; 
        end
        if ~isempty(l)&isnan(x(l(1)))
            gpos=eceftogeod(state.posr(1:3,n));
            if strcmp(prm.trop,'trop_saast')
                x(l(1))=feval(prm.trop,td+t/86400,[0,pi/2],gpos,state.metr(:,n)');
            else
                x(l(1))=0;
            end
            if length(l)>=2, x(l(2))=0; end
            i=[i,l];
        end
        if ~isempty(m)&isnan(x(m(1))) x(m)=0; i=[i,m]; end
    end
end
if ~isempty(i)
    j=find(~isnan(x)); P(i,j)=0; P(j,i)=0; k=sub2ind(size(P),i,i); P(k)=sig(i).^2;
end

% initial satellite orbit parameters -------------------------------------------
function xi=initobtp(n,prm)
srp=prm.sat.srp(n,:)';
switch prm.obt.p_solarpr
case 'srp_simple', xi=[nan;nan;nan;nan;nan;nan]; % Cr*A/mass
case 'srp_rock4',  xi=[1;0;nan;nan;nan;nan];     % [scale,ybias]
case 'srp_gspm',   xi=[1;0;nan;nan;nan;nan];     % [scale,ybias]
case 'srp_gspmm',  xi=[1;0;0;0;;nan;nan];        % [scale,ybias,along,cross]
case 'srp_code',   xi=[srp(4:6);nan;nan;nan];    % [D0,Y0,B0]
case 'srp_code2',  xi=[srp(4:7);nan;nan];        % [D0,Y0,B0,Z0]
case 'srp_code3',  xi=[srp(4:7);0;0];            % [D0,Y0,B0,Z0,X10,X30]
end
if ~isempty(prm.srpprms)
    j=find(strcmp(prm.srpprms(:,1),prm.sats{n}));
    if ~isempty(j), xi=[prm.srpprms{j,2:7}]'; end
end

% initial receiver position ----------------------------------------------------
function xi=initrcvp(n,state,prm)
omge=7.2921151467E-5;
switch prm.rposmodel
case {0,1}, xi=state.posr(1:3,n);
case 2, xi=state.posr(1:6,n);
case 3, dtr=state.clkr(1,n);
        xi=state.U*Rz(omge*dtr)*state.posr(1:3,n);
end

% filter -----------------------------------------------------------------------
function [x,P,phn,vz,res,ds]=filt(td,t,x,ix,P,z,tz,iz,antp,state,phs,pass,prm)
P0=P; vz=zeros(size(z,1),1); res=repmat(nan,size(z,1),3); k=[]; n=1;

% prefit residuals
[v,H,R,phn,i,ds]=measmodel(td,t,x,ix,z,tz,iz,antp,state,phs,prm);
res(i,1)=v; rmsp=rmsx(v);

while ~isempty(v)&n<=prm.maxiter

    % kalman filter update
    [x,P]=filtfun(td,t,x,ix,P0,v,H,R,prm);
    
    % postfit residuals
    [v,H,R,phn,j,ds(i,:)]=measmodel(td,t,x,ix,z(i,:),tz(i),iz(i,:),antp,state,...
                                    phs,prm);
    i=i(j); vz(i)=1; res(i,2:3)=[v,ds(i,12)];
    
    % exclude outliers
    if prm.outlf>0
        if ~isempty(v), k=find(abs(v)>sqrt(diag(R))*prm.outlf); else k=[]; end
    end
    if isempty(k)
        n=n+1;
    else
        vz(i(k))=2; i(k)=[]; v(k)=[]; H(k,:)=[]; R(k,:)=[]; R(:,k)=[]; n=1;
    end
end
if any(diag(P)<0)
    gt_log('filter singlar : t=%s ix=%s',tstr(td,t),sprintf('%d ',find(diag(P)<0)));
end
if ticks(10)
    gmsg('pass-%d-%d-%d : %s n=%4d %4d out=%3d res=%7.4f %7.4f m',pass,...
         tstr(td,t),size(z,1),sum(vz==1),sum(vz==2),rmsp,rmsx(v));
end

% kalman filter update ---------------------------------------------------------
function [x,P]=filtfun(td,t,x,ix,P,v,H,R,prm)
if ~isempty([ix.sata{:},ix.erp,ix.eco])
    i=find(~isnan(x));
    [x(i),P(i,i),vv,stat]=feval(prm.filter,x(i),P(i,i),v,H(:,i),R,prm.outlp);
    if stat
        gt_log('filter error : t=%s err=%d',tstr(td,t),stat);
    end
else
    for n=1:length(prm.rcvs)
        i=ix.rcva{n};
        i=i(~isnan(x(i))); j=find(any(H(:,i),2));
        if ~isempty(j)
            [x(i),P(i,i),vv,stat]=feval(prm.filter,x(i),P(i,i),v(j),H(j,i),...
                                        R(j,j),prm.outlp);
            if stat
                gt_log('filter error : t=%s %s err=%d',tstr(td,t),prm.rcvs{n},stat);
            end
        end
    end
end

% measurement model ------------------------------------------------------------
function [v,H,R,phn,i,ds]=measmodel(td,t,x,ix,z,tz,iz,antp,state,phs,prm)
C=299792458; R=[]; ds=zeros(size(z,1),17);

i=find(~isnan(z(:,1))&(state.ecls(iz(:,1))==0|prm.eclnsigs(iz(:,1))>0));
if size(state.posr,1)<6, state.posr(6,1)=0; end

% ion-free phase measurement model
[y,H,azel,phn,drds,j,ds(i,3:end)]=rangemodel(td,t,x,ix,[tz(i),double(iz(i,:))],...
    state.ephs,state.clks'*C,state.clkr'*C,state.zpdr,state.posr,state.dpos,...
    state.bcpr,state.U,state.dx,state.dy,state.du,state.rsun,prm.ants,antp',...
    phs,prm.elmin,prm.elmax,state.metr,prm.trop,prm.mapf,prm.mpcc,prm.mpcs,...
    prm.f1,prm.f2,prm.corrf,double(strcmp(prm.rcvs,prm.clkref)),...
    double(prm.rattmodel),double(prm.rposmodel==3),state.mfcr);

% residuals (o-c)
i=i(j); v=z(i)-y; tz=tz(i); iz=iz(i,:); ds(i,1:2)=[z(i),y];
if isempty(i), return; end

% measurement error variance
fact=1+(state.ecls(iz(:,1))~=0).*(prm.eclnsigs(iz(:,1))-1);
var=(prm.sig.obs(iz(:,2)).*fact).^2;
if prm.sigweight>0
    var=var.*(1+prm.sigweight*(1./sin(azel(:,2)).^2-1));
end
% add a priori state variences
if prm.sdfact2>0
    for j=1:length(i)
        n=iz(j,1); m=iz(j,2);
        if prm.est.sato(n)>1
            var(j)=var(j)+state.ephv(1:3,n)'*(drds(j,1:3)'*prm.sdfact2).^2;
        end
        if prm.est.satc(n)>1
            var(j)=var(j)+state.clkv(1,n)'*(C*prm.sdfact2)^2;
        end
        if prm.est.rcvp(m)>1
            var(j)=var(j)+state.posv(1:3,m)'*(state.U*drds(j,1:3)'*prm.sdfact2).^2;
        end
    end
end
% differences of observables
if any(strcmp(prm.omodel,{'zerodiff','ppp'}))
    R=diag(var);
else
    [v,H,R,ig]=feval(prm.omodel,v,H,sqrt(var),[tz,double(iz)],prm.ircv);
    ii=[];
    for n=1:length(v), ii(n,1)=i(ig(n,1)==iz(:,1)&ig(n,3)==iz(:,2)); end
    i=ii;
end

% rms --------------------------------------------------------------------------
function y=rmsx(x), if isempty(x), y=0; else y=sqrt(mean(x.^2)); end

% time string ------------------------------------------------------------------
function s=tstr(td,t), s=sprintf('%04d/%02d/%02d %02d:%02d:%02.0f',mjdtocal(td,t));

% check tick -------------------------------------------------------------------
function f=ticks(n)
persistent tp, ts=now; if 86400*(ts-tp)<1/n, f=0; else f=1; tp=ts; end
