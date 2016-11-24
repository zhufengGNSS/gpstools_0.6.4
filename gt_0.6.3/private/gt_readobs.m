function [z,tz,iz,arc,antp,state,abort,stats]=gt_readobs(td,time,state,prm)
%-------------------------------------------------------------------------------
% [system] : GpsTools
% [module] : read clean observation data/navigation messages
% [func]   : read clean observation data/navigation messages
% [argin]  : td,time = day(mjd),time(sec)
%            state   = a priori states
%            prm     = processing parameters (see prm_gpsest_def.m)
% [argout] : z,tz,iz = observation data/observation data time/index
%            arc     = arc info
%                arc(n,1:2) = arc start/end time (sec)
%                arc(n,3:4) = satellite/station index
%                arc(n,5:6) = L1/L2 approx. phase bias (m)
%                arc(n,7:9) = spares
%            antp    = receiver antenna pcv parameters
%                antp(:,:,r) = prm.rcvs{r} antenna pcv
%            state   = a priori states
%                state(n).posr(:,r) = time(n) prm.rcvs{r} guess position/velocity
%                state(n).clkr(1,r) = time(n) prm.rcvs{r} guess clock
%            abort   = abort condition (1:aborted,0:completed)
%            stats   = receiver info. {{rcvtype,anttype,antdel},...}
% [note]   :
% [version]: $Revision: 8 $ $Date: 06/07/13 8:04 $
%            Copyright(c) 2004-2006 by T.Takasu, all rights reserved
% [history]: 06/03/25  0.1  separated from gpsestd.m
%            06/06/24  0.2  add antenna alias
%-------------------------------------------------------------------------------
C=299792458; z=[]; tz=[]; iz=[]; arc=[]; antp=zeros(length(prm.sats),47);
abort=0; stats={}; ch=4+find(strcmp(prm.observ,{'LC','L1','L2'})); % channel

if prm.obs.srfilt
    [nav,inav]=readnav(td,time,prm.sats,prm.rcvs,prm.dirs.nav,prm.src.nav);
    if isempty(nav)
        gt_log('no navigation message   : %s dir=%s',prm.src.nav,prm.dirs.nav);
        return;
    end
end
for n=1:length(prm.rcvs)
    if gmsg('reading clean observation data : %s',prm.rcvs{n}), abort=1; return; end
    
    % read clean observation data
    [zn,izn,q,adel,at,rt,rstat,azel,slip,arcn]=...
        readobs(td,time,prm.sats,prm.rcvs{n},prm.dirs.obc,'obsc',{},inf,prm.tunit,ch);
    
    if ~isempty(zn)
        % exclude satellites/receivers
        tt=round(izn(:,1));
        i=find(time(1)<=tt&tt<=time(end)&mod(tt-time(1),prm.tint)==0);
        
        % exclude satellites/receivers
        for m=1:size(prm.exclude,1)
            t=round(([caltomjd(prm.exclude{m,2}),caltomjd(prm.exclude{m,3})]-td)*86400);
            s=find(strcmp(prm.sats,prm.exclude{m,1}));
            if ~isempty(s)
                i(izn(i,2)==s&t(1)<=tt(i)&tt(i)<=t(2))=[];
            elseif strcmp(prm.rcvs{n},prm.exclude{m,1})
                i(t(1)<=tt(i)&tt(i)<=t(2))=[];
            end
        end
        zn=zn(i,:); izn=izn(i,:); izn(:,3)=n; arcn(:,4)=n;
        
        % extract arcs
        arcn(:,1)=time(1)+ceil((arcn(:,1)-time(1))/prm.tint)*prm.tint;
        arcn(:,2)=time(1)+floor((arcn(:,2)-time(1))/prm.tint)*prm.tint;
        arcn(arcn(:,1)<time(1),1)=time(1);
        arcn(arcn(:,2)>time(end),2)=time(end);
        arcn(arcn(:,1)>=arcn(:,2),:)=[];
        
        % multipath correction with phase residuals
        if prm.obs.srfilt
            zn=correctmp(td,time,zn,izn,prm.rcvs{n},nav,inav,prm);
        end
        z=[z;zn]; tz=[tz;izn(:,1)]; iz=[iz;uint8(izn(:,2:end))];
        arc=[arc;arcn,zeros(size(arcn,1),3)];
        stats{n}={prm.rcvs{n},rt{1},at{1},adel};
        
        % read antenna pcv paremeters
        [apc1,apc2,apv1,apv2,stat]=readpcv(prm.rcvs{n},at{1},prm.pcv);
        if stat
            antp(n,:)=[apc1,apc2,adel([1,3,2])',apv1,apv2];
        else
            gt_log('no rcv antenna pcv      : %s ant=%s file=%s',prm.rcvs{n},at{1},prm.pcv);
        end
        if isempty(zn)
            gt_log('no valid clean obs data : %s dir=%s',prm.rcvs{n},prm.dirs.obc);
        end
    else
        stats{n}={prm.rcvs{n},'','',[0,0,0]};
        gt_log('no clean obs data       : %s dir=%s',prm.rcvs{n},prm.dirs.obc);
    end
    % read receiver guess position/velocity/clock-bias
    rs=repmat(nan,7,length(time));
    if ~isempty(rstat{1})
        [t,i,j]=intersect(time,round(rstat{1}(:,1)));
        rs(:,i)=rstat{1}(j,2:8)';
    end
    for m=1:length(time)
        if prm.est.rcvp(n)<=1, state(m).posr(:,n)=rs(1:6,m); end
        if prm.est.rcvc(n)<=1, state(m).clkr(:,n)=rs(7,m)/C; end
    end
end
if ~isempty(iz)
    [tz,i]=sort(tz); z=z(i,:); iz=iz(i,:);
end

% correct multipath with phase-residuals ---------------------------------------
function z=correctmp(td,time,z,iz,rcv,nav,inav,prm)
sr=86164-9; MU=3.986005E14; corr=zeros(size(z,1),1);
ts=[time(1)+sr*prm.obs.srdays(1)-30,time(end)+sr*prm.obs.srdays(2)+30];
[t,index,ress,outs,sats]=readstats(td,ts,prm.dirs.inpr,'f',prm.tunit,{rcv});
if isempty(t)
    gt_log('no phase-residuals      : %s dir=%s',rcv,prm.dirs.inpr);
    return;
end
for n=1:length(prm.sats)
    i=min(find(inav==n));
    m=find(strcmp(sats,prm.sats{n}));
    if ~isempty(i)&~isempty(m)
        toff=floor(2*2*pi/(sqrt(MU/nav(i,17)^6)+nav(i,12))/prm.tint)*prm.tint;
        i=find(iz(:,2)==n);
        j=find(index(:,1)==m&~outs);
        for k=prm.obs.srdays(1):prm.obs.srdays(2)
            [tt,ii,jj]=intersect(round(iz(i,1)),round(t(j))-k*toff);
            corr(i(ii))=filter([0.5,1,1,0.5],3,ress(j(jj),2));
        end
    end
end
z=z-corr;
