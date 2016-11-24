%-------------------------------------------------------------------------------
% [system] : GpsTools
% [module] : double difference
% [func]   : generate double difference between satellites/stations
% [argin]  : dz  = residuals (O-C) (m)
%            dgds= partial derivatives
%            sig = measurement noise std. deviation (m)
%            ig  = observation index [t,isat,ircv,arcf;...]
%            ircv= station numbers of baselines [ircv1,ircv2;...]
% [argout] : dzs = double diff. of resuduals
%            G   = double diff. of partial derivatives
%            R   = covariance matrix of measurement noise
%            iz  = satellite/station index[isat1,isat2,ircv1,ircv2;...]
% [note]   :
% [version]: $Revision: 12 $ $Date: 2008-11-25 10:02:15 +0900 (火, 25 11 2008) $
%            Copyright(c) 2004-2006 by T.Takasu, all rights reserved
% [history]: 04/06/19  0.1  new
%-------------------------------------------------------------------------------

% (mex function)

