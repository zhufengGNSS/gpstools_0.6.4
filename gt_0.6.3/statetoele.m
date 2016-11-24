%-------------------------------------------------------------------------------
% [system] : GpsTools
% [module] : satellite state to orbit element
% [func]   : convert satellite position/velocity to orbit element
% [argin]  : state = position/velocity(m,m/sec) [x;y;z;xdot;ydot:zdot]
% [argout] : ele   = orbit element(m,deg) [a,e,i,OMG,omg,M]
% [note]   : tangental orbit element
% [version]: $Revision: 2 $ $Date: 06/07/08 14:22 $
%            Copyright(c) 2004-2006 by T.Takasu, all rights reserved
% [history]: 04/01/05   0.1  new
%------------------------------------------------------------------------------

% (mex function)

