function [poss,covs]=readpos(td,time,rcvs,posdir,possrc,tunit,opt)
%-------------------------------------------------------------------------------
% [system] : GpsTools
% [module] : read station positin data
% [func]   : read station coordinates/postions
% [argin]  : td,time = date (mjd-gpst), time vector (sec)
%            rcvs = station list
%           (posdir) = data directory (default:current)
%           (possrc) = data source    (default:'igssnx')
%                   'igssnx'   = IGS Final (weekly solution)
%                   'gsipos'   = GSI position estimation
%                   'itrf2000' = ITRF2000(ITRF2000_GPS.SSC.txt)
%                   'itrf2005' = ITRF2005(ITRF2005_GPS.SSC.txt)
%                   'itrf97'   = ITRF97  (ITRF97_GPS.SSC.txt)
%                   'igs00'    = IGS00   (IGS01P37_RS54.snx)
%                   'igb00'    = IGb00   (IGS03P33_RS99.snx)
%                   'igs05'    = IGS05   (IGS05.snx)
%                   'posf',    = estimated (forward)
%                   'posb'     = estimated (backward)
%                   'posfb'    = estimated (smoothed)
%                   'possf',   = estimated (SINEX,forward)
%                   'possb'    = estimated (SINEX,backward)
%                   'possfb'   = estimated (SINEX,smoothed)
%                   'approx'   = approx. position
%                   ''         = all position data
%           (tunit) = processing unit time (hour)
%           (opts)  = options
%                   'interp'   : interpolate position
% [argout] : poss = coordinate/position (m) [x,y,z] (ecef)
%                  poss(n,:,m) = time(n),rcvs{m} coordinate/position
%            covs = variences (m) [covx,covy,covz]
%                  covs(n,:,m) = time(n),rcvs{m} coordinate/position variences
% [note]   :
% [version]: $Revision: 16 $ $Date: 2008-12-12 15:49:30 +0900 (金, 12 12 2008) $
%            Copyright(c) 2004-2006 by T.Takasu, all rights reserved
% [history]: 04/11/19  0.1  new
%            04/11/23  0.2  use later pos if duplicated value exists
%                           add data source of estimated positions
%            04/11/25  0.3  add geocenter offset correction
%            05/03/04  0.4  add approx. position,igs00,igb00
%                           delete geocenter offset correction
%            05/03/24  0.5  separate readsinexpos.m
%            05/06/13  0.6  separate readgsipos.m
%                           add argin tunit
%            05/07/25  0.7  add possrc:possf,possb,possfb
%            06/03/06  0.8  add argin opts
%            06/06/24  0.9  add warning messages
%            08/11/26  0.10 add itrf2005,igs05 for possrc (gt_0.6.4)
%                           add error message in case of position file absence
%                           fix bug on error-stop if format error
%-------------------------------------------------------------------------------
if nargin<4, posdir=''; end
if nargin<5, possrc=''; end
if nargin<6, tunit=24; end
if nargin<7, opt=''; end
if isempty(possrc), src='igssnx'; else src=possrc; end
if ischar(rcvs), rcvs={rcvs}; end
rcvs=upper(rcvs);

poss=repmat(nan,[length(time),3,length(rcvs)]); covs=poss;
if any(strcmp(possrc,{'posf','posb','posfb'}))
    for m=1:length(rcvs)
        [poss(:,:,m),covs(:,:,m)]=ReadPosEst(td,time,rcvs{m},posdir,possrc,tunit,opt);
    end
elseif any(strcmp(possrc,{'possf','possb','possfb'}))
    for n=1:length(time)
        for m=1:length(rcvs)
            [poss(n,:,m),covs(n,:,m)]=ReadPosEstSinex(td,time(n),rcvs{m},posdir,possrc,tunit);
        end
    end
elseif strcmp(possrc,'gsipos')
    for m=1:length(rcvs)
        for n=1:length(time)
            [poss(n,:,m),covs(n,:,m)]=ReadPosData(td,time(n),rcvs{m},posdir,src,m==1&n==1);
        end
    end
else
    for n=1:length(time)
        for m=1:length(rcvs)
            [poss(n,:,m),covs(n,:,m)]=ReadPosData(td,time(n),rcvs{m},posdir,src,m==1&n==1);
        end
        if isempty(possrc)&any(isnan(poss(n,1,:)))
            for m=find(isnan(poss(n,1,:)))'
                [poss(n,:,m),covs(n,:,m)]=ReadPosData(td,time(n),rcvs{m},posdir,'itrf2000');
            end
        end
        if isempty(possrc)&any(isnan(poss(n,1,:)))
            for m=find(isnan(poss(n,1,:)))'
                [poss(n,:,m),covs(n,:,m)]=ReadPosData(td,time(n),rcvs{m},posdir,'itrf97');
            end
        end
        if isempty(possrc)&any(isnan(poss(n,1,:)))
            for m=find(isnan(poss(n,1,:)))'
                [poss(n,:,m),covs(n,:,m)]=ReadPosData(td,time(n),rcvs{m},posdir,'gsipos');
            end
        end
        if isempty(possrc)&any(isnan(poss(n,1,:)))
            for m=find(isnan(poss(n,1,:)))'
                [poss(n,:,m),covs(n,:,m)]=ReadPosData(td,time(n),rcvs{m},posdir,'approx');
            end
        end
        if isempty(possrc)&any(isnan(poss(n,1,:)))&any(strcmp(possrc,{'igssnx','gsipos'}))
            for m=find(isnan(poss(n,1,:)))'
                gt_log('no position data        : %s dir=%s',rcvs{m},posdir);
            end
        end
    end
end

% read station position data ---------------------------------------------------
function [pos,cov]=ReadPosData(td,ts,rcv,posdir,possrc,reload)
if nargin<6, reload=0; end
persistent td_ rcv_ dir_ src_
persistent epoch time names poss vels psigs vsigs rcvp
if isempty(td_), td_ =0; end
td=td+floor(ts/86400); ts=rem(ts,86400);
pos=[nan;nan;nan]; cov=pos;

switch possrc
case 'igssnx', % IGS Final
    if reload&(td~=td_|~strcmp(posdir,dir_)|~strcmp(possrc,src_))
        gpsd=td+ts/86400-44244; gpsw=floor(gpsd/7);
        epoch=mjdtocal(44244+gpsw*7);
        file=sprintf('igs%02dP%04d.snx',mod(epoch(1),100),gpsw);
        file=gfilepath(posdir,file,epoch);
        [names,poss,vels,psigs,vsigs]=readsinexpos(file);
        if isempty(names)
            gt_log('no position data        : file=%s src=%s',file,possrc);
            return
        end
        td_=td; dir_=posdir; src_=possrc;
    end
    i=find(strcmp(rcv,names)); if isempty(i)|isempty(epoch), return, end
    if length(i)>1
        gt_log('duplicated positions    : %s dir=%s src=%s',rcv,posdir,possrc);
        i=max(i);
    end
    t=(td+ts/86400-caltomjd(epoch))/365.25;
    pos=(poss(:,i)+vels(:,i)*t)';
    cov=(psigs(:,i).^2+vsigs(:,i).^2*t)';

case 'gsipos', % GSI position estimation
    if reload&(td~=td_|~strcmp(rcv,rcv_)|~strcmp(posdir,dir_)|~strcmp(possrc,src_))
        epoch=mjdtocal(td);
        file=sprintf('gsi%s%s.%02d.pos',filesep,rcv,mod(epoch(1),100));
        file=gfilepath(posdir,file,epoch);
        [epoch,time,poss,psigs]=readgsipos(file); poss=poss'; psigs=psigs';
        if isempty(epoch)
            gt_log('no position data        : file=%s src=%s',file,possrc);
            return
        end
        td_=td; rcv_=rcv; dir_=posdir; src_=possrc;
    end
    if isempty(epoch), return, end
    [tdd,tss]=caltomjd(epoch);
    t=(td-tdd)*86400+ts-tss;
    if time(1)<=t&t<=time(end)
        pos=interp1(time,poss',t,'linear','extrap');
    elseif t<time(1)
        pos=poss(:,1)';
    elseif time(end)<t
        pos=poss(:,end)';
    end
    cov=[0,0,0];

case {'itrf2005','itrf2000','itrf97','igs00','igb00','igs05'}, % ITRF,IGS00,IGb00,IGS05
    switch possrc
    case 'itrf2005'
        epoch=[2000,1,1,0,0,0]; type='txt'; file='ITRF2005_GPS.SSC.txt'; 
    case 'itrf2000'
        epoch=[1997,1,1,0,0,0]; type='txt'; file='ITRF2000_GPS.SSC.txt'; 
    case 'itrf97'
        epoch=[1997,1,1,0,0,0]; type='txt'; file='ITRF97_GPS.SSC.txt'; 
    case 'igs00'
        epoch=[1998,1,1,0,0,0]; type='snx'; file='IGS01P37_RS54.snx'; 
    case 'igb00'
        epoch=[1998,1,1,0,0,0]; type='snx'; file='IGS03P33_RS99.snx'; 
    case 'igs05'
        epoch=[2000,1,1,0,0,0]; type='snx'; file='IGS05.snx'; 
    end
    if reload&(~strcmp(posdir,dir_)|~strcmp(possrc,src_))
        file=gfilepath(posdir,file,mjdtocal(td));
        if strcmp(type,'txt')
            [names,poss,vels,psigs,vsigs]=ReadItrfTxt(file);
        else
            [names,poss,vels,psigs,vsigs]=readsinexpos(file);
        end
        if isempty(names)
            gt_log('no position data        : file=%s src=%s',file,possrc);
            return;
        end
        dir_=posdir; src_=possrc;
    end
    i=find(strcmp(rcv,names));
    if ~isempty(i)
        if length(i)>1
            gt_log('duplicated positions    : %s dir=%s src=%s',rcv,posdir,possrc);
            i=max(i);
        end
        t=(td+ts/86400-caltomjd(epoch))/365.25;
        pos=(poss(:,i)+vels(:,i)*t)';
        cov=(psigs(:,i).^2+vsigs(:,i).^2*t)';
    end
case 'approx', % approx. position
    if isempty(rcvp), rcvp=ReadApproxPos; end
    i=find(strcmp(upper(rcv),rcvp(:,1)));
    if ~isempty(i)
        pos=geodtoecef([rcvp{i,[3,2,4]}])';
    end

otherwise,
    gt_log('position data src error : %s src=%s',rcv,possrc); return
end

% read estimated positions -----------------------------------------------------
function [pos,cov]=ReadPosEst(td,time,rcv,posdir,possrc,tunit,opt)
pos=repmat(nan,length(time),3); cov=pos;
[t,xs,vs,prm]=readest(td,time,'pos',rcv,posdir,possrc(4:end),tunit);
if isempty(t), return, end
[tt,i,j]=intersect(time,t);
if isempty(i), return, end
if length(i)==1
    pos=xs(j,:);
    cov=vs(j,:);
elseif strcmp(opt,'interp')
    pos=interp1(tt,xs(j,:),time,'linear','extrap');
    cov=interp1(tt,vs(j,:),time,'linear','extrap');
else
    pos(i,:)=xs(j,:);
    cov(i,:)=vs(j,:);
end

% read estimated positions sinex -----------------------------------------------
function [pos,cov]=ReadPosEstSinex(td,time,rcv,posdir,possrc,tunit)
persistent td_ dir_ src_ names poss psigs
pos=repmat(nan,3,1); cov=pos;
td=td+floor(time/86400);
if isempty(td_)|td~=td_|~strcmp(possrc,src_)|~strcmp(posdir,dir_)
    ep=mjdtocal(td);
    file=sprintf('pos%s%04d%02d%02d.snx',possrc(5:end),ep(1:3));
    file=gfilepath(posdir,file,ep);
    [names,poss,vels,psigs,vsigs]=readsinexpos(file);
    if isempty(names), gt_log('no position data file   : %s',file); return, end
    td_=td; dir_=posdir; src_=possrc;
end
rcv=rcv(max(length(rcv)-3,1):end);
i=find(strcmp(names,rcv)); if isempty(i), return, end
pos=poss(:,i);
cov=psigs(:,i).^2;

% read ITRF text file ----------------------------------------------------------
function [names,poss,vels,psigs,vsigs]=ReadItrfTxt(file)
names={}; poss=[]; vels=[]; psigs=[]; vsigs=[];
fd=fopen(file,'rt'); if fd<0, return, end
if findstr(file,'ITRF2005'), nc=33; else nc=34; end
poss=zeros(3,10000); vels=poss; psigs=poss; vsigs=poss; n=0;
while 1
    str=fgetl(fd); if ~isstr(str), break, end
    if ~isempty(StrToNum(str,1,1))
        name=SubStr(str,nc,4);
        if strcmp(name,'    ')
            if n>0
                for m=1:3
                    vels(m,n)=StrToNum(str,nc-8+m*13,12);
                    vsigs(m,n)=StrToNum(str,nc+38+m*6,5);
                end
            end
        elseif ~isempty(name)
            n=n+1;
            names{n}=name;
            for m=1:3
                poss(m,n)=StrToNum(str,nc-8+m*13,12);
                psigs(m,n)=StrToNum(str,nc+38+m*6,5);
            end
        end
    end
end
poss=poss(:,1:n); vels=vels(:,1:n); psigs=psigs(:,1:n); vsigs=vsigs(:,1:n);
fclose(fd);

% read approx position file ----------------------------------------------------
function rcvp=ReadApproxPos
[dirs,f]=fileparts(which(mfilename));
file=fullfile(dirs,'data','rcvs_poss.txt');
if exist(file)
    [f1,f2,f3,f4]=textread(file,'%s%f%f%f%*[^\n]','delimiter',',',...
                           'commentstyle','matlab');
    rcvp=[f1,num2cell([f2,f3,f4])];
else
    rcvp={};
end

% substring --------------------------------------------------------------------
function sstr=SubStr(str,pos,ns), sstr=str(pos:min(pos+ns-1,length(str)));

% string->number ---------------------------------------------------------------
function num=StrToNum(str,pos,ns)
num=sscanf(str(pos:min(pos+ns-1,length(str))),'%f',1);
if isempty(num), num=nan; end
