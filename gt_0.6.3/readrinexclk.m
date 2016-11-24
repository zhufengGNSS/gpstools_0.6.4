%-------------------------------------------------------------------------------
% [system] : GpsTools
% [module] : read rinex clock data
% [func]   : read rinex clock data file
% [argin]  : file   = file path
% [argout] : epoch  = first epoch [year,month,day,hour,min,sec]
%            time   = time vector relative to epoch
%            types  = clock data types
%            names  = satellite/station list
%            data   = clock data
%                     data(n,:) : index(n) clock data
%            index  = satellite/station index
% [note]   :
% [version]: $Revision: 2 $ $Date: 06/07/08 14:21 $
%            Copyright(c) 2004-2006 by T.Takasu, all rights reserved
% [history]: 04/03/08   0.1  new
%-------------------------------------------------------------------------------

% (mex function)

