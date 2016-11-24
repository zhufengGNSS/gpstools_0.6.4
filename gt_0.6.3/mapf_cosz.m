function [mapfd,mapfw]=mapf_cosz(t,azel,gpos,metprm)
%-------------------------------------------------------------------------------
% [system] : GpsTools
% [module] : tropospheric mapping function - cosz
% [func]   : cosz tropospheric mapping function
% [argin]  : t = date (mjd) (no use)
%            azel  = azimath/elevation angle (rad) [az,el]
%           (gpos) = latitude/longitude/height(deg,m) [lat,lon,h]
%           (metprm) = meteological parameters (no use)
% [argout] : mapfd = dry mapping function
%            mapfw = wet mapping function
% [note]   :
% [version]: $Revision: 3 $ $Date: 06/07/08 14:12 $
%            Copyright(c) 2004-2006 by T.Takasu, all rights reserved
% [history]: 04/09/23  0.1  new
%-------------------------------------------------------------------------------
mapfd=1/sin(azel(2));
mapfw=1/sin(azel(2));
