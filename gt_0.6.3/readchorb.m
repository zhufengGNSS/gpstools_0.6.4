function [epoch,time,ephs]=readchorb(file)
%-------------------------------------------------------------------------------
% [system] : GpsTools
% [module] : read CHORB CHAMP orbit file
% [func]   : read CHORB CHAMP orbit file
% [argin]  : file  = file path
% [argout] : epoch = reference epoch [year,month,day,hour,min,sec]
%            time  = time vector relative to epoch (sec)
%            ephs  = satellite position/velocity (ecef)
%                    ephs(n,1:3) = time(n) satellite position (m)
%                    ephs(n,4:6) = time(n) satellite velocity (m/sec)
% [note]   :
% [version]: $Revision: 2 $ $Date: 06/07/08 14:19 $
%            Copyright(c) 2004-2006 by T.Takasu, all rights reserved
% [history]: 06/04/02  0.1  new
%-------------------------------------------------------------------------------
epoch=[]; time=[]; ephs=[]; n=0;

f=fopen(file,'rt'); if f<0, return, end
while 1
    str=fgetl(f); if ~ischar(str), break, end
    if strcmp(str,'ORBIT'), n=1;
    elseif n>0
        time(n,:)=[str2num(str(1:6))/10,(str2num(str(7:17))-64184000)*1E-6+13];
        for m=1:6, ephs(n,m)=str2num(str(12*m+6:12*m+17)); end
        n=n+1;
    end
end
fclose(f);
if isempty(time), return, end
days=caltomjd([2000,1,1,12,0,0])+time(:,1);
day0=floor(days(1));
epoch=mjdtocal(day0);
time=time(:,2)+(days-day0)*86400;
ephs(:,1:3)=ephs(:,1:3)*1E-3;
ephs(:,4:6)=ephs(:,4:6)*1E-7;
