%-------------------------------------------------------------------------------
% [system] : GpsTools
% [module] : phase windup effect correction
% [func]   : calculate phase windup effect correction
% [argin]  : rsat = satellite position(m) (eci)
%            rsun = sun position(m) (eci)
%            posr = station position (m) (ecef)
%            U    = eci to ecef transformation matrix
%            phwinp = previous correction(rad)
% [argout] : phwin  = phase windup phase correction(rad)
% [note]   :
% [version]: $Revision: 12 $ $Date: 2008-11-25 10:02:15 +0900 (火, 25 11 2008) $
%            Copyright(c) 2004-2006 by T.Takasu, all rights reserved
% [history]: 04/06/18   0.1  new
%-------------------------------------------------------------------------------

% (mex function)

