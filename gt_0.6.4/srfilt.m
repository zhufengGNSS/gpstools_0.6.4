function pos=srfilt(td,time,rcv,pos,srday,srprm,dirs,fb,tunit)
%-------------------------------------------------------------------------------
% [system] : GpsTools
% [module] : sidereal filter
% [func]   : sidereal filter
% [argin]  : td,time = day(mjd), time(sec)
%            rcv     = station
%            pos     = position
%            srday   = relative days/adjustment (day,sec) [start,end,adj]
%            srprm   = bandpass-filter frequencies (hz) [low,high]
%            dirs    = estimated position directory
%            fb      = forward/backward estimation
%            tunit   = estimation unit time (hr)
% [argout] : pos     = corrected position
% [note]   :
% [version]: $Revision: 12 $ $Date: 2008-11-25 10:02:15 +0900 (火, 25 11 2008) $
%            Copyright(c) 2004-2006 by T.Takasu, all rights reserved
% [history]: 06/05/04  0.1  new
%-------------------------------------------------------------------------------
if srday(2)<srday(1), return, end
srsec=86164+srday(3); % sidereal day (sec)
day=srday(1):srday(2);
off=zeros(length(time),3); m=0;
for n=1:length(day)
    t=time+round(srsec*day(n)/(time(2)-time(1)))*(time(2)-time(1));
    corr=readpos(td,t,rcv,dirs,fb,tunit);
    i=find(~isnan(corr(:,1)));
    if ~isempty(i)
        % delete dc
        corr(i,:)=detrend(corr(i,:),'constant');
        
        % interpolate/extrapolate outage period
        if isnan(corr(1,1)), corr(1,:)=corr(i(1),:); i=[1;i]; end
        if isnan(corr(end,1)), corr(end,:)=corr(i(end),:); i=[i;length(t)]; end
        corr=interp1(t(i),corr(i,:),t);
        
        % bandbass filter
        if srprm(1)<srprm(2), corr=bpfilt(t,corr,srprm(1),srprm(2)); end
        
        m=m+1; off=(off*(m-1)+corr)/m;
    end
end
pos=pos-off;
