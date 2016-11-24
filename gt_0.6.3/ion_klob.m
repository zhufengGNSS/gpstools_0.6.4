%-------------------------------------------------------------------------------
% [system] : GpsTools
% [module] : ionospheric delay model - Klobuchar model
% [func]   : calculate ionospheric delay(L1) by Kloubachar model
% [argin]  : t    = date/time(mjd-utc)
%            azel = satellite azimath/elevation angle(rad) [az,el]
%            gpos = latitude/longitude/height(deg,m) [lat,lon,h]
%            (ion_prm) = ionospheric parameters [alpha0,1,2,3;beta0,1,2,3]
%            <global>
%            utc_tai = utc-tai(sec) (default:-32)
% [argout] : ion  = ionospheric delay(L1)(m)
% [note]   : Reference: ICD-GPS-200C Navstar GPS Space Segment/Navigation User
%                       Interfaces Rev.C, 20.3.3.5.2.5
% [version]: $Revision: 2 $ $Date: 06/07/08 14:12 $
%            Copyright(c) 2004-2006 by T.Takasu, all rights reserved
% [history]: 03/12/31   0.1  new
%-------------------------------------------------------------------------------

% (mex function)

