%-------------------------------------------------------------------------------
% [system] : GpsTools
% [module] : geodetic to ecef position
% [func]   : convert geodetic postion to ecef
% [argin]  : gpos = latitude/longitude/height(deg,m) [lat,lon,h] 
% [argout] : epos = ecef position (m) [x;y;z]
%            E    = ecef to local tangental coordinate transformation matrix
% [note]   : WGS84 ellipsoide
% [version]: $Revision: 12 $ $Date: 2008-11-25 10:02:15 +0900 (火, 25 11 2008) $
%            Copyright(c) 2004-2006 by T.Takasu, all rights reserved
% [history]: 04/05/31   0.1  new
%-------------------------------------------------------------------------------

% (mex function)

