%-------------------------------------------------------------------------------
% [system] : GpsTools
% [module] : satellite antenna offset correction
% [func]   : calculate satellite antenna offset correction
% [argin]  : rsat = satellite postion(m) (eci)
%            rrcv = station position(m) (eci)
%            rsun = sun position(m) (eci)
%            apc  = satellite antenna offset(m) [x;y;z] (sallite fixed coordinate)
% [argout] : apcs = satellite antenna offset(m)
% [note]   :
% [version]: $Revision: 12 $ $Date: 2008-11-25 10:02:15 +0900 (火, 25 11 2008) $
%            Copyright(c) 2004-2006 by T.Takasu, all rights reserved
% [history]: 04/06/01   0.1  new
%-------------------------------------------------------------------------------

% (mex function)

