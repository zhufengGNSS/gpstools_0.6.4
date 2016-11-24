%-------------------------------------------------------------------------------
% [system] : GpsTools
% [module] : satellite orbit element to position/velocity
% [func]   : transform satellite orbit elements to state(position/velocity)
%            calculate partial derivatives of state by orbit elements
% [argin]  : ele   = orbit elements(m,deg) [a,e,i,OMG,omg,M]
% [argout] : state = position/velocity(m,m/sec) [x;y;z;xdot;ydot:zdot]
%            (dsde) = partial derivatives of state by orbit elements(6x6)
% [note]   : error if a<=0,e<=0,1<=e
% [version]: $Revision: 2 $ $Date: 06/07/08 14:11 $
%            Copyright(c) 2004-2006 by T.Takasu, all rights reserved
% [history]: 04/01/05   0.1  new
%-------------------------------------------------------------------------------

% (mex function)

