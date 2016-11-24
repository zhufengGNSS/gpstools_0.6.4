function [ephs,sigs]=ephtoecef(td,time,ephs,sigs,prm,sats,dirs,fb,tunit)
%-------------------------------------------------------------------------------
% [system] : GpsTools
% [module] : transform coordinates of satellite ephemeris
% [func]   : transform coordinates of satellite ephemeris (eci to ecef)
% [argin]  : td   = date (mjd-gpst)
%            time = time vector (sec)
%            ephs = satellite ephemeris (eci)
%            sigs = satellite ephemeris standard deviation (eci)
%            prm  = estimation parameter
%            sats = satellites
%            dirs = estimated data directory
%            fb   = forward/backward
%            tunit = processing unit time (hr) (for 'ephf','ephb')
% [argout] : ephs = satellite ephemeris (ecef)
%            sigs = satellite ephemeris standard deviation (ecef)
% [note]   :
% [version]: $Revision: 3 $ $Date: 06/07/08 14:11 $
%            Copyright(c) 2004-2006 by T.Takasu, all rights reserved
% [history]: 05/04/12   0.1  new
%            05/04/29   0.2  support prm.nutmodel
%            05/06/28   0.3  use readerp() for reading erp estimation
%-------------------------------------------------------------------------------
srcs={'','igs','igu','bulb','bula'};
if isfield(prm.src,'erp'), erpsrc=prm.src.erp; else erpsrc=srcs{prm.est.erp}; end
if isfield(prm,'nutmodel'), nutmodel=prm.nutmodel; else nutmodel=''; end

for n=1:length(time)
    tu=td+(time(n)+19+prm.utc_tai)/86400;
    if prm.est.erp==0
        erp_value=prm.erp;
    elseif prm.est.erp==1
        erp_value=readerp(tu,prm.dirs.est,['erp',fb],tunit);
        erp_value(4:5)=prm.erp(4:5);
    else
        erp_value=readerp(tu,prm.dirs.erp,erpsrc);
    end
    if prm.erpvar, erp_value(1:3)=erp_value(1:3)+erpvar(tu,prm.utc_tai); end
    
    [U,a,b,gmst,dx,dy,du]=ecsftoecef(tu,erp_value,prm.utc_tai,nutmodel);
    for m=1:length(sats)
        ephs(n,4:6,m)=ephs(n,4:6,m)*U'+ephs(n,1:3,m)*du';   % vel(ecef)
        ephs(n,1:3,m)=ephs(n,1:3,m)*U';                     % pos(ecef)
        sigs(n,1:3,m)=sqrt(diag(U*diag(sigs(n,1:3,m).^2)*U'))';
        sigs(n,4:6,m)=sqrt(diag(U*diag(sigs(n,4:6,m).^2)*U'))';
    end
end
