%-------------------------------------------------------------------------------
% [system] : GpsTools
% [module] : receiver antenna offset
% [func]   : receiver antenna offset
% [argin]  : azel = azimuth/elevation angle(rad) [az,el]
%            apc1 = L1 phase center offset (m) [up;north;east]
%            apc2 = L2 phase center offset (m) [up;north;east]
%            ecc  = antenna deltas (m) [up;north;east]
%           (apv1)= L1 phase center variation (m) (73x19)
%           (apv2)= L2 phase center variation (m) (73x19)
%                   apv?[i+j*73]=az(i),el(j) pcv (az=0:5:360deg,el=0:5:90deg)
% [argout] : apcr = receiver antenna offset (m) [apcr1,apcr2]
% [note]   :
% [version]: $Revision: 16 $ $Date: 2008-12-12 15:49:30 +0900 (金, 12 12 2008) $
%            Copyright(c) 2004-2008 by T.Takasu, all rights reserved
% [history]: 04/06/01   0.1  new
%            08/11/25   0.4  support pcv azimuth angle dependancy (gt_0.6.4)
%-------------------------------------------------------------------------------

% (mex function)

