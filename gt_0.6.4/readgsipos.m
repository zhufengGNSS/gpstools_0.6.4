%-------------------------------------------------------------------------------
% [system] : GpsTools
% [module] : read GSI position estimation file
% [func]   : read GSI position estimation file
% [argin]  : file = filepath
% [argout] : epoch = epoch time [year,month,day,hour,min,sec]
%            time  = time vector relativ to epoch
%            poss  = position
%                    poss(n,1:3) : time(n) position [x,y,z] (m)
%            psigs = std. deviations
%                    psigs(n,1:3) : time(n) pos std.dev [dx,dy,dz] (m)
% [note]   :
% [version]: $Revision: 16 $ $Date: 2008-12-12 15:49:30 +0900 (金, 12 12 2008) $
%            Copyright(c) 2004-2008 by T.Takasu, all rights reserved
% [history]: 05/06/02  0.1  separated from readpos.m
%-------------------------------------------------------------------------------

% (mex function)

