%-------------------------------------------------------------------------------
% [system] : GpsTools
% [module] : spheric harmonic functions
% [func]   : spheric harmonic functions by azimuth/elevation angle
% [argin]  : az,el = azimuth/elevation angle (rad)
%            nmax  = max degrees of spheric harmonic functions
% [argout] : fc,fs = spheric harmonic functions
%                fc(n+1,m+1) = Pnm(-cos(2*el))*cos(m*az)
%                fs(n+1,m+1) = Pnm(-cos(2*el))*sin(m*az)
%                (Pnm = normalized Legendre polynomial)
% [note]   :
% [version]: $Revision: 12 $ $Date: 2008-11-25 10:02:15 +0900 (火, 25 11 2008) $
%            Copyright(c) 2004-2006 by T.Takasu, all rights reserved
% [history]: 04/10/18  0.1  new
%-------------------------------------------------------------------------------

% (mex function)

