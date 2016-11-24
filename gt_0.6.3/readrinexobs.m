%-------------------------------------------------------------------------------
% [system] : GpsTools
% [module] : read rinex observation data
% [func]   : read rinex observation data file
% [argin]  : file    = file path
% [argout] : epoch   = first obs. epoch [year,month,day,hour,min,sec]
%            time    = time vector relative to epoch (sec)
%            types   = observation data types
%            units   = observation data units
%            sats    = satellites list
%            rcv     = receiving station
%            data    = observation data
%                      data(n,m) : time(n),index(n),types(m) data
%            index   = satellite index
%            rcvpos  = approx station position (m) [x;y;z]
%            antdel  = antenna delta (m) [up;east;north]
%            anttype = antenna model name
%            rcvtype = receiver model name
%            comment = rinex header lines separated by ';'
% [note]   :
% [version]: $Revision: 2 $ $Date: 06/07/08 14:21 $
%            Copyright(c) 2004-2006 by T.Takasu, all rights reserved
% [history]: 04/11/16  0.1  new
%-------------------------------------------------------------------------------

% (mex function)

