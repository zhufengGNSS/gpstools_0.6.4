function [apc,ecc,ant,rec]=readrcv(td,ts,rcv,rcvdir,rcvsrc)
%-------------------------------------------------------------------------------
% [system] : GpsTools
% [module] : read station info
% [func]   : read station infomation
% [argin]  : td,ts = date(mjd-gpst),time(sec)
%            rcv = station
%           (rcvdir) = data directory (default:current)
%           (rcvsrc) = data source    (default:'igssnx')
%                   'igssnx'   = IGS00(sinex)
% [argout] : apc = receiver anntena phase center offset
%                    [up1,up2;north1,north2;east1,east2] (m)
%            ecc = receiver anntena delta [up;north;east] (m)
%            ant = anntena model
%            rec = receiver model
% [note]   :
% [version]: $Revision: 12 $ $Date: 2008-11-25 10:02:15 +0900 (ÁÅ´, 25 11 2008) $
%            Copyright(c) 2004-2006 by T.Takasu, all rights reserved
% [history]: 04/11/19  0.1  new
%-------------------------------------------------------------------------------
persistent file epoch time names eps recs ants apcs eccs, if isempty(file), file=''; end
if nargin<4, rcvdir=''; end
if nargin<5, rcvsrc=''; end
if ~isempty(rcvdir), rcvdir=[rcvdir,filesep]; end
if isempty(rcvsrc), rcvsrc='igssnx'; end
rcv=upper(rcv);

rec=''; ant='';
apc=[0.1100,0.1280;0.0000,0.0000;0.0000,0.0000]; % ÉmÉ~ÉiÉãíl
ecc=[0.0000;0.0000;0.0000];
switch rcvsrc
case 'igssnx', % IGS SINEX
    mjd=td+ts/86400;
    gpsd=mjd-44244; gpsw=floor(gpsd/7);
    dt=MjdToCal(44244+gpsw*7);
    f=sprintf('%sigs%02dP%04d.snx',rcvdir,mod(dt(1),100),gpsw);
    if ~strcmp(f,file)
        fd=fopen(f,'rt');
        if fd<0, disp(['warning : rcv info file open error : ',f]), return, end
        file=f; epoch=dt; time=0;
        [names,eps,recs,ants,apcs,eccs]=ReadSinexRcv(fd);
        fclose(fd);
    end
    if ~isempty(names)
        i=max(find(strcmp(rcv,names)&eps(:,1)<=mjd));
        if ~isempty(i)
            rec=recs{i};
            ant=ants{i};
            apc=[apcs(i,1:3)',apcs(i,4:6)'];
            ecc=eccs(i,:)';
        else
            disp(['warning : rcv info data read error : ',rcv])
        end
    else
        disp(['warning : rcv info data read error : ',rcv])
    end
end

% read receiver/anntena info from sinex file -----------------------------------
function [names,eps,recs,ants,apcs,eccs]=ReadSinexRcv(fd)
eps=zeros(10000,1); apcs=zeros(10000,6); eccs=zeros(10000,3);
start_rec=0; start_ant=0; start_apc=0; start_ecc=0;
while 1
    str=fgetl(fd); if ~isstr(str), break, end

    if findstr(str,'+SITE/RECEIVER'), start_rec=1; n=0;
    elseif findstr(str,'-SITE/RECEIVER'), start_rec=0; end
    
    if findstr(str,'+SITE/ANTENNA'), start_ant=1;
    elseif findstr(str,'-SITE/ANTENNA'), start_ant=0; end
    
    if findstr(str,'+SITE/GPS_PHASE_CENTER'), start_apc=1;
    elseif findstr(str,'-SITE/GPS_PHASE_CENTER'), start_apc=0; end
    
    if findstr(str,'+SITE/ECCENTRICITY'), start_ecc=1;
    elseif findstr(str,'-SITE/ECCENTRICITY'), start_ecc=0; end
    
    if start_rec&SubStr(str,1,1)==' '
        n=n+1;
        eps(n,1)=TimeToMjd(SubStr(str,17,12));
        names{n,1}=SubStr(str,2,4);
        recs{n,1}=deblank(SubStr(str,43,20));
    elseif start_ant&SubStr(str,1,1)==' '
        i=find(strcmp(SubStr(str,2,4),names));
        for ii=i', ants{ii,1}=deblank(SubStr(str,43,20)); end
    elseif start_apc&SubStr(str,1,1)==' '
        ant=deblank(SubStr(str,2,20));
        i=find(strcmp(ant,ants));
        for ii=i', apcs(ii,:)=sscanf(SubStr(str,29,41),'%f')'; end
    elseif start_ecc&SubStr(str,1,1)==' '
        i=find(strcmp(SubStr(str,2,4),names));
        for ii=i', eccs(ii,:)=sscanf(SubStr(str,48,25),'%f')'; end
    end
end
eps=eps(1:n,:); apcs=apcs(1:n,:); eccs=eccs(1:n,:);

% substring --------------------------------------------------------------------
function sstr=SubStr(str,pos,ns), sstr=str(pos:min(pos+ns-1,length(str)));

% sinex time->mjd --------------------------------------------------------------
function mjd=TimeToMjd(str)
mjd=0; t=sscanf(str,'%f:%f:%f');
if length(t)>=3
    if t(1)<=50, t(1)=t(1)+2000; else t(1)=t(1)+1900; end
    mjd=CalToMjd([t(1),1,1])+t(2)-1+t(3)/86400;
end

