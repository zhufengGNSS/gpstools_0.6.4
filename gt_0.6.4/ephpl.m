%-------------------------------------------------------------------------------
% [system] : GpsTools
% [module] : solar/planetary ephemeris
% [func]   : read solar/planetary postions from de405 ephemeris
% [argin]  : t = date/time(mjd-utc)
%            pl= solar/planetary numbers
%                (1:mercury,2:venus,3:earth-moon barycenter,4:mars,5:jupiter,
%                 6:saturn,7:uranus,8:neptune,9:pluto,10:moon,11:sun)
%            <global>
%            utc_tai = utc-tai(sec)
%            ephpdir = ephemris data directory
% [argout] : r = solar/planetary positions(km)
% [note]   : coordinate:barycenter(all except for moon),eci(icrf)(moon)
% [version]: $Revision: 12 $ $Date: 2008-11-25 10:02:15 +0900 (火, 25 11 2008) $
%            Copyright(c) 2004-2006 by T.Takasu, all rights reserved
% [history]: 04/03/13  0.1  new
%-------------------------------------------------------------------------------

% (mex function)

