function [ephs,sigs]=readeph(td,time,sats,ephdir,ephsrc,tunit,opts)
%-------------------------------------------------------------------------------
% [system] : GpsTools
% [module] : read satellite ephemerides
% [func]   : read satellite ephemerides
% [argin]  : td   = date (mjd-gpst)
%            time = time vector (sec)
%            sats = satellites list
%            (ephdir) = input data directory (default:current)
%            (ephsrc) = input data source    (default:'igs')
%                     'igs'  = IGS final
%                     'igr'  = IGS rapid
%                     'igu'  = IGS ultra-rapid
%                     'igp'  = IGS ultra-rapid(predicted half)
%                     'igsc' = IGS final(SP3C)
%                     'igrc' = IGS rapid(SP3C)
%                     'nima','cod','esa',... = NIMA,IGS AC ephemerides
%                     'codr' = code rapid
%                     'ephs' = interpolated by geneph
%                     'ephf' = estimated (forward)
%                     'ephb' = estimated (backward)
%                     'ephfb'= estimated (smoothed)
%                     'brdc' = broadcast ephemeris (IGS Combined)
%            (tunit)  = processing unit time (hr) (for 'ephf','ephb')
%            (opts)   = options
%                     'interp' = interpolate ephemeris data
% [argout] : ephs = satellite positions/velocities
%                ephs(n,:,m)=time(n)sats{m}position/velocy [x,y,z,vx,vy,vz] (m,m/s)
%            sigs = satellite positions/velcities sigma (m)
% [note]   :
% [version]: $Revision: 21 $ $Date: 06/07/08 14:19 $
%            Copyright(c) 2004-2006 by T.Takasu, all rights reserved
% [history]: 04/11/16   0.1  new
%            05/02/28   0.2  delete argin U
%                            separate interpolation to interpeph.m
%            05/04/12   0.3  add 'igp','ephf','ephb','cod',... as ephsrc
%            05/05/20   0.4  add ephsrc 'igsc','igrc'
%                            output sigs in case of sp3 ephemeris
%            05/09/16   0.5  add ephsrc 'brdc'
%            05/09/20   0.6  add argin opts interpolation
%            05/12/08   0.7  improve performance for brdc
%            06/06/24   0.8  add warning messages, add code rapid
%-------------------------------------------------------------------------------
if nargin<4, ephdir=''; end
if nargin<5, ephsrc=''; end
if nargin<6, tunit=24;  end
if nargin<7, opts='';  end
if ischar(sats), sats={sats}; end
if isempty(ephsrc), ephsrc='igs'; end

% read satellite positions (ecef)
switch ephsrc
case {'ephf','ephb','ephfb'}
    [ephs,sigs]=ReadEphEst(td,time,sats,ephdir,ephsrc,tunit);

case 'brdc'
    [ephs,sigs]=ReadEphBrdc(td,time,sats,ephdir);

otherwise
    ephs=repmat(nan,[length(time),6,length(sats)]);
    sigs=zeros([length(time),6,length(sats)]);
    if strcmp(opts,'interp')&~strcmp(ephsrc,'ephs')
        nmax=10; tint=900;
        t=(floor(time(1)/tint)-nmax:floor(time(end)/tint)+nmax)*tint;
        eph=repmat(nan,[length(t),6,length(sats)]); sig=eph;
        for n=1:length(t)
            [eph(n,:,:),sig(n,:,:)]=ReadEphData(td,t(n),sats,ephdir,ephsrc);
        end
        for n=1:length(sats)
            i=find(~isnan(eph(:,1,n)));
            if length(i)>=2
                [ephs(:,1:3,n),ephs(:,4:6,n)]=interplag(t(i),eph(i,1:3,n),time,nmax);
                sigs(:,1:3,n)=interp1(t(i),sig(i,1:3,n),time,'linear','extrap');
                sigs(:,4:6,n)=sigs(:,1:3,n)*sqrt(2)/tint;
            end
        end
    else
        for n=1:length(time)
            [ephs(n,:,:),sigs(n,:,:)]=ReadEphData(td,time(n),sats,ephdir,ephsrc);
        end
    end
end

% read satellite positions -----------------------------------------------------
function [eph,sig]=ReadEphData(td,t,sats,ephdir,ephsrc)
persistent file te ts ti ephs esats sigs, if isempty(file), file=''; end
eph=repmat(nan,6,length(sats)); sig=eph;
switch ephsrc,
case 'ephs'
    dd=floor(t/86400); td=td+dd; t=t-86400*dd; dt=mjdtocal(td,t);
    f=sprintf('%04d%02d%02d%02d',dt(1:3),0);
    d=gfilepath(ephdir,'',dt);
    if ~strcmp([d,f],file)|~strcmp([sats{:}],[esats{:}])
        file=[d,f];
        [epoch,time,ephs,esats,sigs]=ReadEphs(d,f,sats);
        if isempty(epoch), te=0; ts=0; ti=1; return, end
        te=caltomjd(epoch(1:3)); ts=epoch(4)*3600+epoch(5)*60+epoch(6);
        ti=(time(end)-time(1))/(length(time)-1);
    end

case {'igu','igp'} % sp3 ultra-rapid
    if td>=53106, tu=6; else tu=12; end % 6H interval after 2004/4/11
    if strcmp(ephsrc,'igu'), tt=t; else tt=t-3600*tu; ephsrc='igu'; end
    dd=floor(tt/86400); tdd=td+dd; tt=tt-86400*dd; dt=mjdtocal(tdd,tt);
    gpsd=tdd-44244; gpsw=floor(gpsd/7); gpsd=floor(gpsd-gpsw*7);
    fn=sprintf('%s%04d%1d_%02d.sp3',ephsrc,gpsw,gpsd,floor(dt(4)/tu)*tu);
    f=gfilepath(ephdir,fn,dt);
    if ~strcmp(f,file)
        file=f;
        [epoch,time,ephs,esats,accs,sigs]=readsp3(file);
        if all(all(all(isnan(sigs))))
            for n=1:length(esats)
                % if accurate code equal 0, default value used
                if isnan(accs(n)), sigs(:,:,n)=1; else sigs(:,:,n)=accs(n); end
            end
        end
        if size(ephs,2)>=7
            ephs=ephs(:,[1:3,5:7],:);
            sigs=sigs(:,[1:3,5:7],:);
        end
        if isempty(epoch)
            gt_log('no sp3 ephemeris        : %s',file);
            te=0; ts=0; ti=1; return
        end
        te=caltomjd(epoch(1:3)); ts=epoch(4)*3600+epoch(5)*60+epoch(6);
        ti=(time(end)-time(1))/(length(time)-1);
    end

otherwise % sp3 daily
    dd=floor(t/86400); td=td+dd; t=t-86400*dd; dt=mjdtocal(td,t);
    gpsd=td-44244; gpsw=floor(gpsd/7); gpsd=floor(gpsd-gpsw*7);
    switch ephsrc
    case 'nima', fn=sprintf('nim%04d%1d.eph',gpsw,gpsd);
    case 'cod',  fn=sprintf('%s%04d%1d.eph',ephsrc,gpsw,gpsd);
    case 'codr', fn=sprintf('COD%04d%1d.EPH_R',gpsw,gpsd);
    case 'esa',  fn=sprintf('%s%04d%1d.sp3c',ephsrc,gpsw,gpsd);
    case {'igsc','igrc'}, fn=sprintf('%s%04d%1d.sp3c',ephsrc(1:end-1),gpsw,gpsd);
    otherwise,   fn=sprintf('%s%04d%1d.sp3',ephsrc,gpsw,gpsd);
    end
    f=gfilepath(ephdir,fn,dt);
    if ~strcmp(f,file)
        file=f;
        [epoch,time,ephs,esats,accs,sigs]=readsp3(file);
        if all(all(all(isnan(sigs))))
            for n=1:length(esats)
                % if accurate code equal 0, default value used
                if isnan(accs(n)), sigs(:,:,n)=1; else sigs(:,:,n)=accs(n); end
            end
        end
        if size(ephs,2)>=7
            ephs=ephs(:,[1:3,5:7],:);
            sigs=sigs(:,[1:3,5:7],:);
        end
        if isempty(epoch)
            gt_log('no sp3 ephemeris        : %s',file);
            te=0; ts=0; ti=1; return
        end
        te=caltomjd(epoch(1:3)); ts=epoch(4)*3600+epoch(5)*60+epoch(6);
        ti=(time(end)-time(1))/(length(time)-1);
    end
end
for n=1:length(sats)
    j=find(strcmp(sats{n},esats));
    if ~isempty(j)
        i=((td-te)*86400+t-ts)/ti+1;
        if 1<=i&i<=size(ephs,1)&abs(i-floor(i))<1E-3 
            if any(ephs(i,1:3,j)~=0)
                eph(1:3,n)=ephs(i,1:3,j)';
                sig(1:3,n)=sigs(i,1:3,j)';
            end
            if size(ephs,2)>=6
                eph(4:6,n)=ephs(i,4:6,j)';
                sig(4:6,n)=sigs(i,4:6,j)';
            end
        end
    end
end
% read broadcast ephemeris -----------------------------------------------------
function [ephs,sigs]=ReadEphBrdc(td,time,sats,ephdir)
persistent ts ssats nav inav, if isempty(ts), ts=[0,0]; end
if td+time(1)/86400<ts(1)-0.01|ts(2)+1.01<td+time(end)/86400
    ssats={}; for n=1:31, ssats={ssats{:},sprintf('GPS%02d',n)}; end
    [nav,inav]=readnav(td,time,ssats,{},ephdir,'brdc');
    ts=floor(td+[time(1),time(end)]/86400);
end
zoff=zeros(length(sats),1);
%prmsat=prm_gpssats(td);
%[s,i,j]=intersect(sats,prmsat(:,1));
%zoff(i([prmsat{j,2}]==0))=1.023; % zoffset=1.023m for Block II/IIA
%zoff(i([prmsat{j,2}]==1))=1.6;   % zoffset=1.600m for Block IIR

ephs=repmat(nan,[length(time),6,length(sats)]); sigs=ephs;
if ~isempty(nav)
    for n=1:length(sats)
        navs=nav(inav==find(strcmp(sats{n},ssats)),:);
        for m=1:length(time)
            [pos,dts,vel]=navtostate(td,time(m),navs);
            pos=pos+zoff(n)*pos/norm(pos);
            ephs(m,:,n)=[pos',vel'];
        end
    end
else
    gt_log('no navigation messages  : %s dir=%s','brdc',ephdir);
end

% read ephemeris data generated by geneph --------------------------------------
function [epoch,time,ephs,sats,sigs]=ReadEphs(dirs,file,sats)
epoch=[]; time=[]; ephs=[]; sigs=[];
for n=1:length(sats)
    fs=fullfile(dirs,['ephs_',sats{n},'_',file,'.mat']);
    if exist(fs)
        sig=[];
        load(fs)
        ephs(:,:,n)=data;
        if isempty(sig), sigs(:,:,n)=zeros(size(data));
        else sigs(:,:,n)=sig; end
    else
        if ~isempty(time), ephs(:,:,n)=nan; sigs(:,:,n)=nan; end
        gt_log('no interpolated ephems  : %s',fs);
    end
end

% read estimated ephemeris in processing unit time -----------------------------
function [ephs,sigs]=ReadEphEst(td,time,sats,dirs,src,tunit)
ephs=repmat(nan,[length(time),6,length(sats)]); sigs=ephs;
for n=1:length(sats)
    [t,xs,vs,prm]=readest(td,time,'eph',sats{n},dirs,src(4:end),tunit);
    [t,i,j]=intersect(time,t);
    if ~isempty(i)
        ephs(i,:,n)=xs(j,1:6);
        sigs(i,:,n)=sqrt(vs(j,1:6));
    end
end
if ~isempty(prm)
    [ephs,sigs]=ephtoecef(td,time,ephs,sigs,prm,sats,dirs,src(4:end),tunit);
end

% coordinate rotation matrix ---------------------------------------------------
function R=Rz(t), R=[cos(t),sin(t),0;-sin(t),cos(t),0;0,0,1];
