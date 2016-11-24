%-------------------------------------------------------------------------------
% [system] : GpsTools
% [module] : satellite-station range
% [func]   : calculate range between satellite and station
% [argin]  : state = satellite position/velocity (m) (eci)
%            rpos  = station position (m) (ecef)
%            U     = eci to ecef transformation matrix
%            dtr   = receiver clock bias (sec)
%            (corrlightt) = light-time/receiver clock correction flag (1:on)
%            (direc) = ranging direction (1:sat to sta,2:sta to sat)
% [argout] : rs    = satellite tx position(m) (eci)
%            range = range between satellite and station(m)
%            (drds) = partical derivatives of range by satellite state
% [note]   :
% [version]: $Revision: 2 $ $Date: 06/07/14 16:16 $
%            Copyright(c) 2004-2006 by T.Takasu, all rights reserved
% [history]: 04/04/12   0.1  new
%-------------------------------------------------------------------------------

% (mex function)

