%-------------------------------------------------------------------------------
% [system] : GpsTools
% [module] : navigation messages to satellite position/velocity
% [func]   : calculate satllite position/velocity from navigation messages
% [argin]  : td    = date (mjd-gpst)
%            ts    = time (sec)
%            nav   = navigation messages
%           (type) = GPS/GLONASS('N'=GPS,'G'=GLONASS)
% [argout] : pos = satellite position [x;y;z] (m)(ecef)
%            dts = satellite clock error [bias;drift;drift-rate]
%                  bias/drift/drift-rate(sec,sec/sec,sec/sec^2)
%           (vel)= satellite velocity [vx;vy;vz] (m/sec)(ecef)
% [note]   : no TGD correction
% [version]: $Revision: 2 $ $Date: 06/07/08 14:12 $
%            Copyright(c) 2004-2006 by T.Takasu, all rights reserved
% [history]: 04/11/16  0.1  new
%-------------------------------------------------------------------------------

% (mex function)

