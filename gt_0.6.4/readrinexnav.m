%-------------------------------------------------------------------------------
% [system] : GpsTools
% [module] : read rinex navigation message
% [func]   : read rinex navigation message file
% [argin]  : file    = file path
% [argout] : sats    = satellites list
%            rcv     = receiving station
%            data    = navigation messages
%                      data(n,:) : index(n) navigation message
%            index   = satellite index
%            ionprm  = ionospheric parameters
%                      ionprm(1,:) : ion alpha
%                      ionprm(2,:) : ion beta
%            dutc    = delta utc parameters [A0,A1,T,W]
%            comment = rinex header separated by ';'
% [note]   :
% [version]: $Revision: 12 $ $Date: 2008-11-25 10:02:15 +0900 (火, 25 11 2008) $
%            Copyright(c) 2004-2006 by T.Takasu, all rights reserved
% [history]: 04/11/17  0.1  new
%-------------------------------------------------------------------------------

% (mex function)

