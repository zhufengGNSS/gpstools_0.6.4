function [dzs,G,R,iz]=zerodiff(dz,dgds,sig,ig,ircvs)
%-------------------------------------------------------------------------------
% [system] : GpsTools
% [module] : zero difference (undifference)
% [func]   : generate undifferenced measurement
% [argin]  : dz  = resuduals
%            dgds= partial derivatives of model
%            sig = meas. noise std. deviations
%            ig  = meas. index [t,isat,ircv,arcf;...]
%            ircvs = baselines of stations [ircv1,ircv2,...] (no use)
% [argout] : dzs = resuduals
%            G   = partial derivatives
%            R   = meas. noise covariences
%            iz  = meas index [isat,0,ircv,0;...]
% [note]   :
% [version]: $Revision: 2 $ $Date: 06/07/08 1:16 $
%            Copyright(c) 2004-2006 by T.Takasu, all rights reserved
% [history]: 04/11/02  0.1  new
%-------------------------------------------------------------------------------
dzs=dz;
G=dgds;
R=diag(sig.^2);
iz=ig(:,[2,2,3,3]); iz(:,[2,4])=0;
