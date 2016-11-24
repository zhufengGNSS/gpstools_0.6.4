%-------------------------------------------------------------------------------
% [system] : GpsTools
% [module] : tropospheric model - GPT
% [func]   : global pressure and temperature based spherical harmonics
% [argin]  : t    = date/time (mjd)
%            gpos = latitude/longitude/height(deg,m) [lat,lon,h]
% [argout] : pres = pressure (hPa)
%            temp = temperature (C)
%            undu = geoid undulation (m)
% [note]   : reference :
%            J. Boehm, R. Heinkelmann, and H. Schuh (2007), Short Note: A global
%            model of pressure and temperature for geodetic applications,
%            Journal of Geodesy , doi:10.1007/s00190-007-0135-3
%            gpt.f (fortran source code) :
%            http://mars.hg.tuwien.ac.at/~ecmwf1/
% [version]: $Revision: 16 $
%            Copyright(c) 2004-2008 by T.Takasu, all rights reserved
% [history]: 08/11/24   0.1  new
%-------------------------------------------------------------------------------

% (mex function)
