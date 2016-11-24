%-------------------------------------------------------------------------------
% [system] : GpsTools
% [module] : eci to satellite-fixed coordinate transformation matrix
% [func]   : calculate eci to satellite-fixed coordinate transformation matrix
% [argin]  : state = satellite position/velocity[x;y;z;vx;vy;vz](m) (eci)
% [argout] : E = eci to satellite-fixed coordinate transformation matrix(3x3)
%                [r_radial;r_along-track;r_cross-track]=E*r
% [note]   :
% [version]: $Revision: 12 $ $Date: 2008-11-25 10:02:15 +0900 (火, 25 11 2008) $
%            Copyright(c) 2004-2006 by T.Takasu, all rights reserved
% [history]: 04/02/08   0.1  new
%-------------------------------------------------------------------------------

% (mex function)

