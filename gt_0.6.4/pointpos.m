%--------------------------------------------------------------------------------
% [system] : GpsTools
% [module] : point positioning
% [argin]  : td   = date(mjd-gpst)
%            ts   = time(sec)
%            z    = observation data(pseudo-range)
%            iz   = observation data satellite index
%            nav  = navigation messages
%            inav = navigation messages satellite/station index
%            posr = approx station position (m) (ecef)
%            cdtr = approx receiver clock bias (m)
% [argout] : posr = station postion (m) (ecef)
%            cdtr = receiver clock bias (m)
% [note]   :
% [version]: $Revision: 12 $ $Date: 2008-11-25 10:02:15 +0900 (火, 25 11 2008) $
%            Copyright(c) 2004-2006 by T.Takasu, all rights reserved
% [history]: 04/10/19  0.1  new
%-------------------------------------------------------------------------------

% (mex function)

