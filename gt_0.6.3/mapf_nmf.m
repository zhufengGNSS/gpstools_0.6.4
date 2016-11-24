%-------------------------------------------------------------------------------
% [system] : GpsTools
% [module] : tropospheric mapping function - NMF
% [func]   : calculate tropospheric mapping function by Niell mapping function
% [argin]  : t    = date/time (mjd)
%            azel = azimath/elevation angle(rad) [az,el]
%            gpos = latitude/longitude/height(deg,m) [lat,lon,h]
%            (metprm) = meteological parameters(no use)
% [argout] : mapfd = dry mapping function
%            mapfw = wet mapping function
% [note]   : reference :
%            A.E.Niell, Global mapping functions for the atmosphere delay at
%            radio wavelengths, 1996
% [version]: $Revision: 2 $ $Date: 06/07/08 14:12 $
%            Copyright(c) 2004-2006 by T.Takasu, all rights reserved
% [history]: 04/05/31   0.1  new
%-------------------------------------------------------------------------------

% (mex function)

