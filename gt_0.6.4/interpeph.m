function ephs=interpeph(td,t,poss,sats,time,erpsrc,dirs)
%-------------------------------------------------------------------------------
% [system] : GpsTools
% [module] : interpolate satellite ephemeris
% [func]   : interpolate satellite ephemeris
% [argin]  : td   = date (mjd-gpst)
%            t    = satellite position time (sec)
%            poss = satellite position (m) (ecef)
%            sats = satellite lists
%            time = interpolated time vector (sec)
%            erpsrc = erp source
%            dirs = data directories
% [argout] : ephs = interpolated ephemeris (m,m/sec) (ecef)
% [note]   :
% [version]: $Revision: 12 $ $Date: 2008-11-25 10:02:15 +0900 (火, 25 11 2008) $
%            Copyright(c) 2004-2006 by T.Takasu, all rights reserved
% [history]: 05/02/27   0.1  separated from readeph.m
%            05/03/07   0.2  fix bug on time systems
%            05/06/26   0.3  delete argin utc_tai, add argin erpsrc
%-------------------------------------------------------------------------------
if ischar(sats), sats={sats}; end

% transform coordinates (ecef->eci)
for n=1:length(t)
    utc_tai=prm_utc_tai(td+t(n)/86400,1);
    tutc=td+(t(n)+19+utc_tai)/86400;
    erp=readerp(tutc,dirs.erp,erpsrc);
    U=ecsftoecef(tutc,erp,utc_tai);
    for m=1:length(sats)
        poss(n,:,m)=poss(n,:,m)*U; % ecef->eci
    end
end
% interpolate satellite position and velocity (eci)
vels=repmat(nan,length(sats),3);
for n=1:length(t)-1
    i=find(t(n)<=time&time<=t(n+1));
    for m=1:length(sats)
        [ephs(i,:,m),vels(m,:)]=...
            interpe(td,t(n:n+1),poss(n:n+1,:,m),vels(m,:),time(i),sats{m});
    end
end
% transform coordinates (eci->ecef)
for n=1:length(time)
    utc_tai=prm_utc_tai(td+time(n)/86400,1);
    tutc=td+(time(n)+19+utc_tai)/86400;
    erp=readerp(tutc,dirs.erp,erpsrc);
    [U,P,N,g,dx,dy,du]=ecsftoecef(tutc,erp,utc_tai);
    for m=1:length(sats)
        ephs(n,4:6,m)=ephs(n,4:6,m)*U'+ephs(n,1:3,m)*du';
        ephs(n,1:3,m)=ephs(n,1:3,m)*U';
    end
end

% interpolate positions and velocities -----------------------------------------
function [ephs,vel]=interpe(td,t,pos,vel,time,sat);
ephs=repmat(nan,length(time)+2,6);
ephs([1,end],1:3)=pos;
if any(isnan(pos(:,1))), ephs=ephs(2:end-1,:); return, end
if isnan(vel(1)), ephs(1,4:6)=(pos(2,:)-pos(1,:))/(t(2)-t(1));
else ephs(1,4:6)=vel; end
for n=1:10
    [ephs,phi]=predorbit(td,[t(1);time(:);t(2)],ephs(1,:),sat);
    dpos=ephs(end,1:3)-pos(2,:); if norm(dpos)<1E-3, break, end 
    ephs(1,4:6)=(ephs(1,4:6)'-phi(1:3,:)\dpos')';
end
if n==10, disp(['warning : interpolation divergence : ',sat]); end
ephs=ephs(2:end-1,:);
vel=ephs(end,4:6);

% propagate satellite orbit ----------------------------------------------------
function [ephs,phi]=predorbit(td,t,eph0,sat)
ephs=zeros(length(t),6); ephs(1,:)=eph0;
for n=1:length(t)-1
    if t(n+1)==t(n)
        ephs(n+1,:)=ephs(n,:);
    else
        utc_tai=prm_utc_tai(td+t(n)/86400,1);
        tutc=td+(t(n)+19+utc_tai)/86400;
        ephs(n+1,:)=state_precorbit(tutc,t(n+1)-t(n),ephs(n,:)',sat)';
    end
end
dv=diag([0.001,0.001,0.001]);
utc_tai=prm_utc_tai(td+t(1)/86400,1);
tutc=td+(t(1)+19+utc_tai)/86400;
for n=1:3
    ephe=state_precorbit(tutc,t(end)-t(1),ephs(1,:)'+[0;0;0;dv(:,n)],sat);
    phi(:,n)=(ephe-ephs(end,:)')/dv(n,n); % phi=dx/dv0
end
