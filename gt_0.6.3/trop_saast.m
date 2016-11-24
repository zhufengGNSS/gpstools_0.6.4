%-------------------------------------------------------------------------------
% [system] : GpsTools
% [module] : tropospheric model - Saastamoinen model
% [func]   : calculate tropospheric delay by Saastamoinen model
% [argin]  : t    = date/time (mjd-utc)(no use)
%            azel = azimuth/elevation angle(rad) [az,el]
%            gpos = latitude/longitude/height(deg,m) [lat,lon,h]
%           (met_prm) = meteorological parameters
%              met_prm(1) = pressure(hPa)
%              met_prm(2) = temperture(C)
%              met_prm(3) = relative humidity(%)
%              (default:standard atmosphere,humidity:50%)
% [argout] : tropd = total/dry tropospheric delay (m)
%           (tropw)= wet tropospheric delay (m)
% [note]   : Reference: H.Lichtenegger GPS Theory and Practice 6.3
% [version]: $Revision: 2 $ $Date: 06/07/08 14:22 $
%            Copyright(c) 2004-2006 by T.Takasu, all rights reserved
% [history]: 03/12/12   0.1  new
%-------------------------------------------------------------------------------

% (mex function)

