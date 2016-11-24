%-------------------------------------------------------------------------------
% [system] : GpsTools
% [module] : navigation messages to satellite position/velocity
% [func]   : calculate satllite position/velocity from navigation messages
% [argin]  : td    = date (mjd-gpst)
%            ts    = time (sec)
%            nav   = navigation messages
%           (type) = GPS/GLONASS('N'=GPS,'G'=GLONASS)
%           (opt)  = option (1:relativity correction off)
% [argout] : pos = satellite position [x;y;z] (m)(ecef)
%            dts = satellite clock error [bias;drift;drift-rate]
%                  bias/drift/drift-rate(sec,sec/sec,sec/sec^2)
%           (vel)= satellite velocity [vx;vy;vz] (m/sec)(ecef)
%           (svh)= sv health
% [note]   : no TGD correction
% [version]: $Revision: 16 $ $Date: 2008-12-12 15:49:30 +0900 (金, 12 12 2008) $
%            Copyright(c) 2004-2008 by T.Takasu, all rights reserved
% [history]: 04/11/16  0.1  new
%            06/02/28  0.2  add argin opt
%            08/12/04  0.3  add argout svh (gt_0.6.4)
%-------------------------------------------------------------------------------

% (mex function)

