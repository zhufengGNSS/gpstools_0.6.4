function [epoch,time,poss,psigs]=readgsipos(file)
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
% [version]: $Revision: 2 $ $Date: 06/07/08 14:21 $
%            Copyright(c) 2004-2006 by T.Takasu, all rights reserved
% [history]: 05/06/02  0.1  separated from readpos.m
%-------------------------------------------------------------------------------
epoch=[]; time=[]; poss=[]; psigs=[];
fd=fopen(file,'rt'); if fd<0, return, end
time=zeros(366,1); poss=zeros(366,1); psigs=poss; start=0; n=0;
while 1
    str=fgetl(fd); if ~isstr(str), break, end
    if findstr(str,'+DATA'), start=1;
    elseif findstr(str,'-DATA'), break; end
    if start
        t=sscanf(SubStr(str,2,19),'%d %d %d %d:%d:%f',6);
        if length(t)==6
            if isempty(epoch), epoch=t'; end
            n=n+1;
            time(n)=(caltomjd(t')-caltomjd(epoch))*86400;
            for m=1:3, poss(n,m)=StrToNum(str,4+18*m,17); end
        end
    end
end
time=time(1:n); poss=poss(1:n,:); psigs=psigs(1:n,:);
fclose(fd);

% substring --------------------------------------------------------------------
function sstr=SubStr(str,pos,ns), sstr=str(pos:min(pos+ns-1,length(str)));

% string->number ---------------------------------------------------------------
function num=StrToNum(str,pos,ns)
num=sscanf(str(pos:min(pos+ns-1,length(str))),'%f');
if isempty(num), num=nan; end
