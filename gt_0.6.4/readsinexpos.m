function [rcvs,poss,vels,psigs,vsigs,goff,gsig]=readsinexpos(file)
%-------------------------------------------------------------------------------
% [system] : GpsTools
% [module] : read sinex station position file
% [func]   : read sinex station position file
% [argin]  : file = sinex station position file
% [argout] : rcvs = station names
%            poss = position (m) [x,y,z;...]
%            vels = velocity (m/year) [vx,vy,vz;...]
%            psigs= position std dev. (m) [dx,dy,dz;...]
%            vsigs= velocity std dev. (m/year) [dvx,dvy,dvz;...]
%            goff = geocenter offset (m) [x,y,z]
%            gsig = geocenter offset std dev. (m) [dx,dy,dz;...]
% [note]   :
% [version]: $Revision: 16 $ $Date: 2008-12-12 15:49:30 +0900 (金, 12 12 2008) $
%            Copyright(c) 2004-2008 by T.Takasu, all rights reserved
% [history]: 05/03/24  0.1  separated from readpos.m
%            08/11/26  0.2  fix bug on error-stop if format error
%-------------------------------------------------------------------------------
rcvs={}; poss=[]; vels=[]; psigs=[]; vsigs=[];
fd=fopen(file,'rt'); if fd<0, return, end
poss=zeros(3,2000); vels=poss; psigs=poss; vsigs=poss; section=0; n=0;
goff=[0;0;0]; gsig=[0;0;0];
while 1
    str=fgetl(fd); if ~isstr(str), break, end
    if     findstr(str,'+SOLUTION/ESTIMATE'), section=1;
    elseif findstr(str,'-SOLUTION/ESTIMATE'), break; end
    if section==1
        label=SubStr(str,8,4);
        switch label,
        case 'STAX',
            n=n+1; rcvs{n}=SubStr(str,15,4);
            poss(1,n)=StrToNum(str,48,21); psigs(1,n)=StrToNum(str,70,11);
        case 'STAY',
            if strcmp(rcvs{n},SubStr(str,15,4))
                poss(2,n)=StrToNum(str,48,21); psigs(2,n)=StrToNum(str,70,11);
            end
        case 'STAZ',
            if strcmp(rcvs{n},SubStr(str,15,4))
                poss(3,n)=StrToNum(str,48,21); psigs(3,n)=StrToNum(str,70,11);
            end
        case 'VELX',
            if strcmp(rcvs{n},SubStr(str,15,4))
                vels(1,n)=StrToNum(str,48,21); pvels(1,n)=StrToNum(str,70,11);
            end
        case 'VELY',
            if strcmp(rcvs{n},SubStr(str,15,4))
                vels(2,n)=StrToNum(str,48,21); pvels(2,n)=StrToNum(str,70,11);
            end
        case 'VELZ',
            if strcmp(rcvs{n},SubStr(str,15,4))
                vels(3,n)=StrToNum(str,48,21); pvels(3,n)=StrToNum(str,70,11);
            end
        case 'XGC ', goff(1)=StrToNum(str,48,21); gsig(1)=StrToNum(str,70,11);
        case 'YGC ', goff(2)=StrToNum(str,48,21); gsig(2)=StrToNum(str,70,11);
        case 'ZGC ', goff(3)=StrToNum(str,48,21); gsig(3)=StrToNum(str,70,11);
        end
    end
end
poss=poss(:,1:n); vels=vels(:,1:n); psigs=psigs(:,1:n); vsigs=vsigs(:,1:n);
fclose(fd);

% substring --------------------------------------------------------------------
function sstr=SubStr(str,pos,ns), sstr=str(pos:min(pos+ns-1,length(str)));

% string->number ---------------------------------------------------------------
function num=StrToNum(str,pos,ns)
num=sscanf(str(pos:min(pos+ns-1,length(str))),'%f',1);
if isempty(num), num=nan; end
