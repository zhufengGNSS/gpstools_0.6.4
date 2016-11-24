%-------------------------------------------------------------------------------
% [system] : GpsTools
% [module] : station displacement 
% [func]   : calculate station displacement by earth tides
% [argin]  : tut     = date/time(mjd-utc)
%            posr    = station position (m) (ecef)
%            possun  = sun position (m) (ecef)
%            posmoon = moon position (m) (ecef)
%            odisp   = ocean tide amplitudes (m)
%            ophas   = ocean tide phases (deg)
%            gmst    = Greenwich mean sidereal time (rad)
%            xpyp    = pole offset(xp,yp) (rad)
%            opt     = option flags [solid,oload,polar,perm] (1:on)
% [argout] : dpos    = station displacement (m) (ecef)
% [note]   :
% [version]: $Revision: 12 $ $Date: 2008-11-25 10:02:15 +0900 (火, 25 11 2008) $
%            Copyright(c) 2004-2006 by T.Takasu, all rights reserved
% [history]: 04/10/19  0.1  new
%-------------------------------------------------------------------------------

% (mex function)

