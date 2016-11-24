%-------------------------------------------------------------------------------
% [system] : GpsTools
% [module] : eci to ecef transformation matrix
% [func]   : calculate eci to ecef transformation matrix
% [argin]  : tm = date/time (mjd-utc)
%            (erp_value)= earth rotation parameters (default:[0,0,0,0,0])
%                erpvalue(1:2) = pole offset xp,yp(rad)
%                erpvalue(3)   = ut1-utc(sec)
%                erpvalue(4:5) = nutation offset dpsi,deps(rad)
%            (utc_tai)  = utc-tai(sec) (default:-32)
%            (model_nut)= nutation model(no use)
% [argout] : U    = eci to ecef transformation matrix(3x3)
%            P,N  = precession,nutation matrix(3x3)
%            gmst = Greenwich mean sidereal time(rad)
%            (dUdxp,dUdyp,dUddt) = dU/dxp,dU/dyp,dU/d(UT1-UTC)(1/rad,1/rad,1/sec)
% [note]   : iau1976 precession,iau1980 nutation model
%            eci  : icrf (J2000.0 mean equator coordinate)
%            ecef : itrf
%            inv(U)=U',inv(P)=P',inv(N)=N'
% [version]: $Revision: 12 $ $Date: 2008-11-25 10:02:15 +0900 (火, 25 11 2008) $
%            Copyright(c) 2004-2006 by T.Takasu, all rights reserved
% [history]: 03/12/31  0.1  new
%-------------------------------------------------------------------------------

% (mex function)

