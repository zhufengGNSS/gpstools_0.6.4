function [clks,sigs,ref]=readclk(td,time,names,clkdir,clksrc,tunit,varargin)
%-------------------------------------------------------------------------------
% [system] : GpsTools
% [module] : read clock data
% [func]   : read clock data
% [argin]  : td,time = date (mjd-gpst),time vector (sec)
%            names = satellites/stations list
%           (clkdir) = clock data directory (default:current)
%           (clksrc) = clock data source (default:'igs')
%                    'igs'  = IGS final (rinex clk)
%                    'igr'  = IGS rapid (rinex clk)
%                    'igu'  = IGS ultra-rapid (sp3)
%                    'igss' = IGS final (sp3)
%                    'igrs' = IGS rapid (sp3)
%                    'cod','emr','esa',... = IGS AC
%                    'codr' = CODE rapid
%                    'clkf' = estimated clock (forward)
%                    'clkb' = estimated clock (backward)
%                    'clkfb'= estimated clock (smoothed)
%                    'brdc' = broadcast clock (IGS Combined)
%                    'igscod'  = IGS final/CODE combined
%                    'igrcodr' = IGS Rapid/CODE repid combined
%                    'igs30s'  = IGS final 30s (rinex clk)
%                    'cod5s'   = CODE 5s (rinex clk)
%                    'igscod5s'= IGS final/CODE 5s combined
%           (tunit) = processing unit time (hr) (for 'clkf','clkb','clkfb')
%           (opts)  = options
%                    'interp'       = interpolate clock data
% [argout] : clks = clock data
%                clks(n,:,m)=time(n),names{m} clock bias [dt] (sec)
%            sigs = clock data standard deviation
%            ref  = reference clock
% [note]   :
% [version]: $Revision: 16 $ $Date: 2008-12-12 15:49:30 +0900 (金, 12 12 2008) $
%            Copyright(c) 2004-2008 by T.Takasu, all rights reserved
% [history]: 04/11/02  0.1  new
%            05/05/31  0.2  return std. dev. in case of clksrc='igs','igr'
%            06/02/27  0.3  add 'brdc','igscod' as clksrc, add argin tunit
%            06/04/22  0.4  add IGS/estimated combined clock
%            06/06/24  0.5  add warning messages, add code rapid
%            07/03/26  0.6  fix bug to read igs ultra-rapid (gt_0.6.3p5)
%            07/03/30  0.7  add 'igs30s' as clksrc (gt_0.6.3p5)
%            08/11/21  0.8  fix bug of error on reading broadcast ephemeris (gt_0.6.4)
%                           support cod5s, igscod5s as clksrc
%-------------------------------------------------------------------------------
if nargin<4, clkdir=''; end
if nargin<5, clksrc=''; end
if nargin<6, tunit=24; end
if nargin<7, opts={}; end
if ischar(names), names={names}; end
if isempty(clksrc), clksrc='igs'; end
interp=0; n=1; 
while n<=nargin-6
    switch varargin{n}
    case 'interp'
        if any(strcmp(clksrc,{'igu','igp','igss','igrs'})), tint=900; else tint=300; end
        time=[time(:);(floor(time(end)/tint)+1)*tint]; interp=1;
    end
    n=n+1;
end
names=upper(names); ref='';
clks=repmat(nan,[length(time),1,length(names)]);
sigs=zeros([length(time),1,length(names)]);
switch clksrc
case 'brdc'
    [clks,sigs]=ReadClkBrdc(td,time,names,clkdir);

case {'clkf','clkb','clkfb'}
    for n=1:length(names)
        [clks(:,:,n),sigs(:,:,n),ref]=ReadEstClk(td,time,names{n},clkdir,clksrc,tunit);
    end
case 'igscod'
    [clks,sigs]=ReadClkCombined(td,time,names,clkdir,{'igs','cod'});

case 'igrcodr'
    [clks,sigs]=ReadClkCombined(td,time,names,clkdir,{'igr','codr'});

case 'igscod5s'
    [clks,sigs]=ReadClkCombined(td,time,names,clkdir,{'igs','cod5s'});

otherwise
    ts=time(1);
    while ts<=time(end)
        [t,ts,clk,sig,name]=ReadClkData(td,ts,clkdir,clksrc);
        [tt,i,j]=intersect(time,t);
        [nn,k,m]=intersect(names,name);
        clks(i,1,k)=clk(j,m); sigs(i,1,k)=sig(j,m);
    end
end
if interp
    for n=1:length(names)
        i=find(~isnan(clks(:,1,n)));
        if length(i)>=2&length(i)<length(time)
            clks(:,1,n)=interp1(time(i),clks(i,1,n),time,'linear');
            sigs(:,1,n)=interp1(time(i),sigs(i,1,n),time,'linear');
            for j=find(diff(time(i))>900)'
                k=find(time(i(j))<time&time<time(i(j+1)));
                clks(k,1,n)=nan; sigs(k,1,n)=nan;
            end
        end
    end
    clks=clks(1:end-1,:,:);
    sigs=sigs(1:end-1,:,:);
end

% read clock data --------------------------------------------------------------
function [time,tn,clk,sig,name]=ReadClkData(td,ts,clkdir,clksrc)
time=[]; clk=[]; sig=[]; name={}; tn=(floor(ts/86400)+1)*86400;
switch clksrc
case {'igu','igp'} % sp3 ultra-rapid
    if td>=53106, tu=6; else tu=12; end % 6H interval after 2004/4/11
    if strcmp(clksrc,'igu'), tt=ts; else tt=ts-3600*tu; clksrc='igu'; end
    dd=floor(tt/86400); tdd=td+dd; tt=tt-86400*dd; dt=mjdtocal(tdd,tt);
    gpsd=tdd-44244; gpsw=floor(gpsd/7); gpsd=floor(gpsd-gpsw*7);
    file=sprintf('%s%04d%1d_%02d.sp3',clksrc,gpsw,gpsd,floor(dt(4)/tu)*tu);
    file=gfilepath(clkdir,file,dt);
    [epoch,time,ephs,name]=readsp3(file);
    tn=(floor(ts/tu/3600)+1)*tu*3600;
    if isempty(epoch)
        gt_log('no sp3 ephemeris/clock  : %s',file);
        return;
    end
    tdd=caltomjd(epoch(1:3)); t=epoch(4:6)*[3600;60;1];
    time=time+t+(tdd-td)*86400;
    clk=squeeze(ephs(:,4,:)); sig=zeros(size(clk));

case {'igss','igrs'} % sp3 daily (final,rapid)
    dd=floor(ts/86400); tdd=td+dd; ts=ts-86400*dd; dt=mjdtocal(tdd,0);
    gpsd=tdd-44244; gpsw=floor(gpsd/7); gpsd=floor(gpsd-gpsw*7);
    file=sprintf('%s%04d%1d.sp3',clksrc(1:end-1),gpsw,gpsd);
    file=gfilepath(clkdir,file,dt);
    [epoch,time,ephs,name]=readsp3(file);
    if isempty(epoch)
        gt_log('no sp3 ephemeris/clock  : %s',file);
        return;
    end
    tdd=caltomjd(epoch(1:3)); t=epoch(4:6)*[3600;60;1];
    time=time+t+(tdd-td)*86400;
    clk=squeeze(ephs(:,4,:)); sig=zeros(size(clk));

otherwise % rinex clk daily
    dd=floor(ts/86400); tdd=td+dd; ts=ts-86400*dd; dt=mjdtocal(tdd,0);
    gpsd=tdd-44244; gpsw=floor(gpsd/7); gpsd=floor(gpsd-gpsw*7);
    switch clksrc
    case 'igs30s', file=sprintf('igs%04d%1d.clk_30s',gpsw,gpsd);
    case 'codr', file=sprintf('COD%04d%1d.CLK_R',gpsw,gpsd);
    case 'cod5s',file=sprintf('COD%04d%1d.CLK_05S',gpsw,gpsd);
    otherwise,   file=sprintf('%s%04d%1d.clk',clksrc,gpsw,gpsd);
    end
    file=gfilepath(clkdir,file,dt);
    [epoch,time,type,name,data,index,sigs]=readrinexclk(file);
    if isempty(epoch)
        gt_log('no rinex clock          : %s',file);
        return;
    end
    tdd=caltomjd(epoch(1:3)); t=epoch(4:6)*[3600;60;1];
    [time,i,j]=unique(time+t+(tdd-td)*86400);
    sigs(isnan(sigs))=0;
    clk=repmat(nan,length(time),length(name));
    sig=zeros(length(time),length(name));
    i=sub2ind(size(clk),j,index);
    clk(i)=data; sig(i)=sigs;
end

% read combined clock ----------------------------------------------------------
function [clks,sigs]=ReadClkCombined(td,time,names,clkdir,srcs)
clks=repmat(nan,[length(time),1,length(names)]); sigs=clks;
ts=[floor(time(1)/300),ceil(time(end)/300)]*300;
if strcmp(srcs{2},'cod5s'), ti=5; else ti=30; end
t1=ts(1):300:ts(2); t2=ts(1):ti:ts(2);
[clks1,sigs1]=readclk(td,t1,names,clkdir,srcs{1}); % interval=300s
[clks2,sigs2]=readclk(td,t2,names,clkdir,srcs{2}); % interval=30s
if isempty(clks1)|isempty(clks2), return, end
warning off
for n=1:length(names)
    [t,i,j]=intersect(t2,t1);
    if length(i)>1
        offset=clks2(:,1,n)-interp1(t2(i),clks2(i,1,n),t2)';
        k=find(~isnan(offset));
        if ~isempty(k)
            offset(abs(offset)>5*std(offset(k)))=nan; % exclude outlier (5-sigma)
        end
        clks3=interp1(t1(j),clks1(j,1,n),t2)'+offset;
        [t,i,j]=intersect(time,t2);
        clks(i,1,n)=clks3(j);
        sigs(:,1,n)=interp1(t1,sigs1(:,1,n),time);
    else
        clks(i,1,n)=clks1(j,1,n);
        sigs(i,1,n)=sigs1(j,1,n);
    end
end
warning on

% read estimated clock data ----------------------------------------------------
function [clk,sig,ref]=ReadEstClk(td,time,name,clkdir,clksrc,tunit)
C=299792458; ref='';
clk=repmat(nan,length(time),1); sig=clk;
[t,xs,vs,prm]=readest(td,time,'clk',name,clkdir,clksrc(4:end),tunit);
[tt,i,j]=intersect(time,t);
if ~isempty(i)
    clk(i,1)=xs(j,1)/C;
    k=find(vs(j,1)>=0); sig(i(k),1)=sqrt(vs(j(k),1))/C;
    if ~isempty(prm), ref=prm.clkref; end
end

% read broadcast ephemeris -----------------------------------------------------
function [clks,sigs]=ReadClkBrdc(td,time,sats,clkdir)
persistent ts ssats nav inav, if isempty(ts), ts=[0,0]; end
if td+time(1)/86400<ts(1)-0.01|ts(2)+1.01<td+time(end)/86400
    ssats={}; for n=1:64, ssats={ssats{:},sprintf('GPS%02d',n)}; end
    [nav,inav]=readnav(td,time,ssats,{},clkdir,'brdc');
    ts=floor(td+[time(1),time(end)]/86400);
end
clks=repmat(nan,[length(time),1,length(sats)]); sigs=clks;
if ~isempty(nav)
    for n=1:length(sats)
        i=find(strcmp(sats{n},ssats));
        if ~isempty(i)
            navs=nav(inav==i,:);
            for m=1:length(time)
                [pos,dts,vel]=navtostate(td,time(m),navs,'N',1);
                clks(m,1,n)=dts(1);
            end
        end
    end
else
    gt_log('no navigation messages  : %s dir=%s','brdc',clkdir);
end
