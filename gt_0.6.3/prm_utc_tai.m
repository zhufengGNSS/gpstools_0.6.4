function utc_tai=prm_utc_tai(t,tsys)
%-------------------------------------------------------------------------------
% [system] : GpsTools
% [module] : UTC-TAI parameters
% [func]   : UTC-TAI parameters
% [argin]  : (t)    = time in mjd (default:recent value)
%            (tsys) = time system (0:UTC,1:GPST) (default:UTC)
% [argout] : utc_tai = UTC-TAI (sec)
% [note]   : valid after 1972/1/1
% [version]: $Revision: 5 $ $Date: 06/07/08 14:17 $
%            Copyright(c) 2004-2006 by T.Takasu, all rights reserved
% [history]: 05/06/09  0.1  new
%            05/09/07  0.2  add leap second on 2005/12/31
%            06/07/07  0.3  separate utc_tai parameter file
%-------------------------------------------------------------------------------
persistent utc_tai_hist

if nargin<1, t=99999999; end
if nargin<2, tsys=0; end

if isempty(utc_tai_hist)
    [dirs,f]=fileparts(which(mfilename));
    file=fullfile(dirs,'data','utc_tai.txt');
    if exist(file)
        [f1,f2]=textread(file,'%f%f','delimiter',',','commentstyle','matlab');
        utc_tai_hist=[f1,f2];
    else
        utc_tai_hist=[-33,0]; % default value (2006/1/1-)
    end
end
if tsys
    i=find(utc_tai_hist(:,2)-utc_tai_hist(:,1)/86400<=t+19/86400); % t in GPST
else
    i=find(utc_tai_hist(:,2)<=t); % t in UTC
end
if ~isempty(i), utc_tai=utc_tai_hist(max(i),1); else utc_tai=0; end
