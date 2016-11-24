function h=plottrop(varargin)
%-------------------------------------------------------------------------------
% [system] : GpsTools
% [module] : plot tropospheric parameters
% [func]   : plot tropospheric parameters
% [argin]  : 'prm',prm    : plot parameter struct
% [argout] : h = figure handle
% [note]   :
% [version]: $Revision: 23 $ $Date: 06/07/13 18:52 $
%            Copyright(c) 2004-2006 by T.Takasu, all rights reserved
% [history]: 04/11/22  0.1  new
%-------------------------------------------------------------------------------
if nargin>0&strncmp(varargin{1},'cb',2), feval(varargin{:})
else h=PlotMain(varargin{:}); end

% main -------------------------------------------------------------------------
function h=PlotMain(varargin)
h=gut('newfigg','','Tropospheric Parameters',[600,400,0,-68],[mfilename,' cbClose']);
if isempty(h), return, end
prm=loadprm('prm_plottrop','prm_plottrop_def');
if nargin>=2&strcmp(varargin{1},'prm'), prm=varargin{2}; end
cb1=[mfilename,' cbPlot1'];
cb2=[mfilename,' cbPlot2'];
cb3=[mfilename,' cbPlot3'];
gut('newmenu','data','&Data ',...
    {'&Read...','&New Window','-','Out&put...','&Export Plot...','-',...
     '&Map Area','&Options...'},{},...
    {[mfilename,' cbRead'],[mfilename,' cbNew'],'',[mfilename,' cbOut'],...
     [mfilename,' cbExport'],'',[mfilename,' cbMap'],[mfilename,' cbOpt']});
gut('newmenu','plot','&Plot',...
    {'&Zenith Total Delay','&Gradient East/North','&ZTD Error','-',...
     'ZTD Error by &Stations','ZTD Error by &Dates','-','ZTD Map','-','&Summary...'},...
    {'ztd','trg','ztde','','errsat','errdate','','map','',''},...
    {cb1,cb1,cb1,'',cb2,cb2,'',cb2,'',cb3});
gut('newmenu','rcvs','&Station',prm.rcvs,prm.rcvs,[mfilename,' cbRcv']);
gut('chkmenu','plot',prm.type);
gut('chkmenu','rcvs',prm.rcv);
set(gut('newbtnh','sbtn',[0,0,36,16],{'<','>'},[mfilename,' cbSwBtn']),...
    'fontsize',8);
set(gcf,'resizefcn',[mfilename,' cbResize']);
LocatePos;
data.prm=prm; data.trop={}; set(h,'userdata',data);
if length(varargin)>0, ReadData; UpdatePlot; end

% callback on new --------------------------------------------------------------
function cbNew
data=get(gcf,'userdata');
feval(mfilename,'prm',data.prm);

% callback on close ------------------------------------------------------------
function cbClose
SaveSetting;
closereq

% callback on resize -----------------------------------------------------------
function cbResize, LocatePos;

% callback on menu read --------------------------------------------------------
function cbRead
sel1={{'Estimated(Forward)','Estimated(Backward)','Estimated(Smoothed)'},...
      {'zpdf','zpdb','zpdfb'}};
sel2={{'IGS Final','IGS Monthly','JMA MSM Online','JMA GSM Online','Saastamoinen'},...
      {'igssnx','igsmon','mso','gso','saast'}};
data=get(gcf,'userdata'); prm=data.prm;
prm1={
't','Start Time (GPST)',[]
't','End Time (GPST)',  []
'n','Time Interval (sec)', 300
'n','Processing Unit Time (hr)',  24
'b','Stations',[mfilename,' cbSelRcv']
'p','Tropos Parameters',{{sel1{1}{:},sel2{1}{:}},{sel1{2}{:},sel2{2}{:}}}
'p','Reference Tropos Parameters',{{'',sel2{1}{:}},{'',sel2{2}{:}}}
};
prm2={
' ','Tropospheric Parameters Directory',''
'd','',prm.dirs.est
' ','Reference Tropos Parameters Directory',''
'd','',prm.dirs.trop
};
gut('newdlg','','Read Data',[314,294]);
gut('newprms','prm1',[12,114,295,25,118],prm1);
gut('newprms','prm2',[12,34,295,18],prm2);
gut('newokcancelbtn','',[127,4,180,23]);
gut('setprms','prm1',prm.tstart,prm.tend,prm.tint,prm.tunit,'',prm.fb,prm.ref);
gut('setudata','prm1_5',prm.rcvs);
if ~gut('waitok'), return, end
[prm.tstart,prm.tend,prm.tint,prm.tunit,q,prm.fb,prm.ref]=gut('getprms','prm1');
[q,prm.dirs.est,q,prm.dirs.trop]=gut('getprms','prm2');
prm.rcvs=gut('getudata','prm1_5');
if ~isempty(prm.rcvs), prm.rcv=prm.rcvs{1}; end
close
data.prm=prm; set(gcf,'userdata',data)
SaveSetting; ReadData; UpdatePlot;

% callback on menu select station ---------------------------------------------
function cbSelRcv
[rcvs,ok]=editlist('Stations',get(gcbo,'userdata'),'rcv');
if ok, set(gcbo,'userdata',rcvs); end

% callback on menu output results ----------------------------------------------
function cbOut
data=get(gcf','userdata'); prm=data.prm;
i=find(strcmp(prm.rcv,prm.rcvs));
if isempty(data.trop)|isempty(i), return, end
s1={};
s1=adds(s1,'%%  TROPOSPHERIC PARAMETER : %s : %s-%s GPST',prm.rcvs{i},tstr(prm.td,prm.time(1)),tstr(prm.td,prm.time(end)));
s1=adds(s1,'%%  date     time(sec)    ztd(m)   grade(m)  gradn(m)   sdztd(m)  sdge(m)   sdgn(m)');
s1=adds(s1,'%%---------------------------------------------------------------------------------');
s2=cell(1,length(prm.time));
for n=1:length(prm.time)
    d=mjdtocal(prm.td,prm.time(n)); t=mod(prm.time(n),86400);
    s2{n}=sprintf('%04d,%02d,%02d,%9.2f, %9.5f,%9.5f,%9.5f, %9.5f,%9.5f,%9.5f\n',...
                  d(1:3),t,data.trop{i}(n,1:3),data.sig{i}(n,1:3));
end
Viewer({s1{:},s2{:}},['Satellite Orbit : ',prm.dirs.est,' (',prm.fb,')']);

% callback on menu export plot -------------------------------------------------
function cbExport
persistent file dirs
if isempty(file), file='plottrop'; end, if isempty(dirs), dirs=''; end
f=gcf; data=get(f,'userdata'); prm=data.prm;
sel1={{'EPS (*.eps)','JPEG (*.jpg)','TIFF (*.tif)','META FILE (*.emf)'},...
      {{'-depsc','-tiff'},{'-djpeg','-r0'},{'-dtiff','-r0'},{'-dmeta'}}};
prm1={
'e','File Name',file; 'p','File Format',sel1; ' ','Directory','';'d','',dirs
};
gut('newdlg','','Export Plots',[274,140]);
gut('newprms','prm1',[12,34,255,25,140],prm1);
gut('newokcancelbtn','',[88,4,180,23]);
if ~gut('waitok'), return, end
[file,format,q,dirs]=gut('getprms','prm1');
close
fs=fullfile(dirs,file); opts={'-noui'};
if ~isempty(strmatch('zt',prm.type))
    rcv=prm.rcv;
    for n=1:length(prm.rcvs)
        prm.rcv=prm.rcvs{n}; data.prm=prm; set(f,'userdata',data)
        figure(f), UpdatePlot;
        print(format{:},opts{:},sprintf('%s%02d',fs,n));
    end
    prm.rcv=rcv; data.prm=prm; set(f,'userdata',data); UpdatePlot;
else
    print(format{:},opts{:},fs);
end

% callback on menu options -----------------------------------------------------
function cbOpt
sel1={{'OFF','ON'},0:1};
sel2={{'Dot','Line','Line and Dot'},1:3};
data=get(gcf,'userdata'); prm=data.prm;
prm1={
'p','Show Reference Parameters',       sel1
'p','Show Mean and RMS Errors',        sel1
'n','Show Est. Std Dev. (sigma,0:OFF)', 0
'p','Plot Style',                      sel2
's','Plot Line Width / Marker Size',   [0.5,5]
's','Plot Axis Range (m)',             [-0.5,0.5]
};
q='';
gut('newdlg','','Options',[320,178]);
gut('newprms','prm1',[12,34,300,23,120],prm1);
gut('newokcancelbtn','',[133,4,180,23]);
gut('setprms','prm1',prm.showf(1),prm.showf(2),prm.nsig,prm.ptype,prm.psize,prm.range)
if ~gut('waitok'), return, end
[prm.showf(1),prm.showf(2),prm.nsig,prm.ptype,prm.psize,range]=gut('getprms','prm1');
if range(1)<range(2), prm.range=range; else prm.range=[-0.5,0.5]; end
close
data.prm=prm; set(gcf,'userdata',data)
SaveSetting; UpdatePlot;

% callback on menu map area ----------------------------------------------------
function cbMap
data=get(gcf,'userdata'); prm=data.prm;
[prm.map,ok]=editmap('Map Area',prm.map);
if ok, data.prm=prm; set(gcf,'userdata',data); UpdatePlot; end

% callback on menu plot --------------------------------------------------------
function cbPlot1
types={'ztd','ztg','ztde'};
data=get(gcf,'userdata'); prm=data.prm;
gut('togglechk',gcbo);
gut('unchkmenu','plot',{'errsat','errdate','map'});
type=gut('getchk','plot');
prm.type={};
for n=1:length(type)
    if any(strcmp(type{n},types))
        if strcmp(type{n},'ztg')
            if size(data.trop{1},2)>=3
                prm.type={prm.type{:},'ztge','ztgn'};
            elseif size(data.trop{1},2)>=6
                prm.type={prm.type{:},'ztgee','ztgen','ztgnn'};
            end
        else
            prm.type={prm.type{:},type{n}};
        end
    end
end
data.prm=prm; set(gcf,'userdata',data)
UpdatePlot

function cbPlot2
data=get(gcf,'userdata'); prm=data.prm;
gut('unchkmenu','plot');
prm.type=get(gcbo,'userdata');
gut('chkmenu','plot',prm.type);
data.prm=prm; set(gcf,'userdata',data)
UpdatePlot

function cbPlot3
data=get(gcf,'userdata'); prm=data.prm;
prm.sum=get(gcbo,'userdata');
data.prm=prm; set(gcf,'userdata',data)
UpdateSummary

% callback on menu station -----------------------------------------------------
function cbRcv
data=get(gcf,'userdata'); prm=data.prm;
prm.rcv=get(gcbo,'userdata');
gut('unchkmenu','rcvs');
gut('chkmenu','rcvs',prm.rcv);
data.prm=prm; set(gcf,'userdata',data);
UpdatePlot

% callback on push swithch button ----------------------------------------------
function cbSwBtn
data=get(gcf,'userdata'); prm=data.prm;
j=find(strcmp(prm.rcvs,prm.rcv));
switch get(gcbo,'tag')
case 'sbtn_1', prm.rcv=prm.rcvs{max(j-1,1)};
case 'sbtn_2', prm.rcv=prm.rcvs{min(j+1,length(prm.rcvs))};
end
gut('unchkmenu','rcvs'); gut('chkmenu','rcvs',prm.rcv);
data.prm=prm; set(gcf,'userdata',data);
UpdatePlot

% locate position --------------------------------------------------------------
function LocatePos
p=get(gcf,'position'); m=15;
gut('setpos','sbtn_1',[p(3)-m-32,p(4)-16]);
gut('setpos','sbtn_2',[p(3)-m-16,p(4)-16]);

% save settingsrs --------------------------------------------------------------
function SaveSetting
data=get(gcf,'userdata'); prm=data.prm;
[path,file]=fileparts(which(mfilename));
save(fullfile(path,'settings','prm_plottrop.mat'),'prm');

% update plot ------------------------------------------------------------------
function UpdatePlot
data=get(gcf,'userdata'); prm=data.prm;
if any(strcmp(prm.type,'errsat')), PlotEstErrRcv
elseif any(strcmp(prm.type,'errdate')), PlotEstErrDate
elseif any(strcmp(prm.type,'map')), PlotTropMap
else PlotTrop; end

% plot tropos parameters -------------------------------------------------------
function PlotTrop
data=get(gcf,'userdata'); prm=data.prm;
if isempty(data.trop), return, end
i=find(strcmp(prm.rcv,prm.rcvs));
trop=data.trop{i}; sig=data.sig{i}; err=data.err{i}; ref=data.ref{i};

ti=['Tropospheric Parameters ',prm.rcv,' : ',tstr(prm.td,prm.time(1)),'-',tstr(prm.td,prm.time(end))];
clf
margin=[0.08,0.03,0.04,0.012];
for n=1:length(prm.type)
    if n<length(prm.type), topts='nolabel'; else topts=''; end
    ggt('subplotv','est',n,length(prm.type),'move',[1,1,1,1],...
        'link',[1,0,1,0],'taxis',prm.td,'topts',topts,'margin',margin,...
        'xlim',[prm.time(1),prm.time(end)]/3600);
    
    switch prm.type{n}
    case 'ztd'
        if prm.showf(1)
            h=plot(prm.time/3600,ref,'m.','markersize',prm.psize(2));
            if ~isempty(prm.ref), legend(h,RefName(prm.ref)); end
        end
        plotest(prm.time,trop(:,1),sig(:,1),prm.nsig,prm.ptype,prm.psize);
        range=meann(trop(:,1))+prm.range;
        if isinf(range(1)), range(1)=max(err(:,1)); end
        if isinf(range(2)), range(2)=min(err(:,1)); end
        if any(isnan(range)), range=[0,1]; end
        set(gca,'ylim',range);
        ylabel('ZTD (m)')
    case 'ztge'
        if size(trop,2)>=2
            plotest(prm.time,trop(:,2),sig(:,2),prm.nsig,prm.ptype,prm.psize);
            ylabel('Grad East (m)')
        end
    case 'ztgn'
        if size(trop,2)>=3
            plotest(prm.time,trop(:,3),sig(:,3),prm.nsig,prm.ptype,prm.psize);
            ylabel('Grad North (m)')
        end
    case 'ztgee'
        if size(trop,2)>=4
            plotest(prm.time,trop(:,4),sig(:,4),prm.nsig,prm.ptype,prm.psize);
            ylabel('Grad EE (m)')
        end
    case 'ztgen'
        if size(trop,2)>=5
            plotest(prm.time,trop(:,5),sig(:,5),prm.nsig,prm.ptype,prm.psize);
            ylabel('Grad EN (m)')
        end
    case 'ztgnn'
        if size(trop,2)>=6
            plotest(prm.time,trop(:,6),sig(:,6),prm.nsig,prm.ptype,prm.psize);
            ylabel('Grad NN (m)')
        end
    case 'ztde'
        plotest(prm.time,err(:,1),[],0,prm.ptype,prm.psize);
        if prm.nsig>0
            plot(prm.time/3600,-prm.nsig*sig(:,1),'r:')
            plot(prm.time/3600, prm.nsig*sig(:,1),'r:')
        end
        ylabel('ZTD Error (m)')
        set(gca,'ylim',prm.range);
        if prm.showf(2)
            s=sprintf('REF: %s  MEAN: %6.4fm RMS: %6.4fm',RefName(prm.ref),...
                      meann(err),rmsn(err));
            ggt('mtext',['est_',num2str(n)],s,1);
        end
    end
    if n==1, title(ti); end
end

% plot tropos parameters in map -----------------------------------------------
function PlotTropMap
data=get(gcf,'userdata'); prm=data.prm;
if isempty(data.trop), return, end
i=find(strcmp(prm.rcv,prm.rcvs));
trop=data.trop{i}; sig=data.sig{i}; err=data.err{i}; ref=data.ref{i};
[fn,fs]=gut('getfont','g');
pos=[0.01,0.012,0.9,0.92];
clf, axis off, box on
gmt('mmap','proj',prm.map.proj,'cent',prm.map.cent,'base',prm.map.base,...
    'scale',prm.map.scale,'pos',pos,'fontname',fn,'fontsize',fs);
gmt('mcoast','lcolor',prm.map.color{1},'scolor',prm.map.color{2},'ccolor',prm.map.color{3});
gmt('mgrid','gint',prm.map.gint(1),'lint',prm.map.gint(2),'color',prm.map.color{4});
[xl,yl]=gmt('getlim');
range=(prm.range(2)-prm.range(1))/10;
label=[num2str(range),'m'];
ti=[tstr(prm.td,prm.time(1)),'-',tstr(prm.td,prm.time(end))];
ti=['Tropospheric Parameters : ',ti];
cs=[2.2,2.7];
cm=get(gcf,'colormap'); nc=size(cm,1); s=(cs(end)-cs(1))/(nc-1);
for n=1:length(prm.rcvs)
    gpos=eceftogeod(readpos(prm.td,0,prm.rcvs{n},'','approx')');
    m=floor((data.trop{n}(1,1)-cs(1))/s)+1; if m<1, m=1; elseif m>nc, m=nc; end
    if ~isnan(m), gmt('mplot',gpos(2),gpos(1),cm(m,:),'marker','.','markersize',prm.psize(2)); end
end
title(ti);
ggt('colorbarv',[0.93,0.025,0.015,0.9],cs,'','fontname',fn,'fontsize',fs);

% plot estimation error by stations -------------------------------------------
function PlotEstErrRcv
data=get(gcf,'userdata'); prm=data.prm;
clf
pos=[0.08,0.12,0.89,0.81];
for n=1:length(data.err) err(n,:)=rmsn(data.err{n}(:,1)); end
err=[err;meann(err)]; range=[0,prm.range(2)];
h=ggt('barplot',err,{prm.rcvs{:},'Average'},'ylim',range,'position',pos);
legend(h,'ZTD')
ylabel('ZTD RMS Error (m)');
ti=['Tropospheric Parameters : ',tstr(prm.td,prm.time(1)),'-',tstr(prm.td,prm.time(end)),...
    ' (REF: ',RefName(prm.ref),')'];
title(ti);

% plot estimatin error by date -------------------------------------------------
function PlotEstErrDate
data=get(gcf,'userdata'); prm=data.prm;
clf
pos=[0.08,0.12,0.89,0.81];
tu=prm.tunit*3600;
ts=(floor(prm.time(1)/tu):floor(prm.time(end)/tu))*tu;
for n=1:length(ts)
    i=find(ts(n)<=prm.time&prm.time<=ts(n)+tu-prm.tint); e=[];
    for m=1:length(data.err)
        if ~isempty(data.err{m}), e=[e;rmsn(data.err{m}(i,1))]; end
    end
    err(n,:)=meann(e);
    dt=mjdtocal(prm.td,ts(n));
    labels{n}=sprintf('%d/%d',dt(2:3));
end
err=[err;meann(err)]; range=[0,prm.range(2)];
h=ggt('barplot',err,{labels{:},'Average'},'ylim',range,'position',pos);
legend(h,'ZTD')
ylabel('ZTD RMS Error (m)');
ti=['Troposphric Parameters : ',tstr(prm.td,prm.time(1)),'-',tstr(prm.td,prm.time(end)),...
    ' (REF: ',RefName(prm.ref),')'];
title(ti);

% update summary ---------------------------------------------------------------
function UpdateSummary
data=get(gcf','userdata'); prm=data.prm;
s={};
s=adds(s,'              TROPOSPHERIC PARAMETERS : ESTIMATED ZTD ERROR');
s=adds(s,' TIME      : %s-%s GPST',tstr(prm.td,prm.time(1)),tstr(prm.td,prm.time(end)));
s=adds(s,' REFERENCE : %s',RefName(prm.ref));
s=AddSummary(s,data,prm,prm.time(1),prm.time(end));
s=adds(s,'');
s=adds(s,' -----------------------------------------------------------------------------');
tu=prm.tunit*3600;
for ts=(floor(prm.time(1)/tu):floor(prm.time(end)/tu))*tu
    s=adds(s,'');
    s=AddSummary(s,data,prm,ts,ts+tu-prm.tint);
end
if isempty(prm.files), ti=prm.dirs.est; else ti=prm.files{1}; end
Viewer(s,ti);

% add summary string -----------------------------------------------------------
function s=AddSummary(s,data,prm,ts,te)
data=get(gcf','userdata'); prm=data.prm;
s=adds(s,' TIME      : %s-%s GPST',tstr(prm.td,ts),tstr(prm.td,te));
s=adds(s,'               mean(m)    rms(m)     # of est    # of ref ');
s=adds(s,' -----------------------------------------------------------------------------');
p=[]; v=[]; r=[];
i=find(ts<=prm.time&prm.time<=te);
for n=1:length(prm.rcvs)
    p(n,1)=meann(data.err{n}(i,1)); v(n,1)=rmsn(data.err{n}(i,1));
    r(n,1)=length(find(~isnan(data.trop{n}(i,1))));
    rr(n,1)=length(find(~isnan(data.ref{n}(i,1))));
    s=adds(s,' %-7s : %10.4f %10.4f %10.0f  %10.0f',prm.rcvs{n},p(n),v(n),r(n),rr(n));
end
s=adds(s,' -----------------------------------------------------------------------------');
s=adds(s,' %-7s : %10.4f %10.4f %10.0f  %10.0f','average',meann(p),meann(v),meann(r),meann(rr));
s=adds(s,' %-7s : %10.4f %10.4f %10.0f  %10.0f','median',mediann(p),mediann(v),mediann(r),meann(rr));

% read tropos data -------------------------------------------------------------
function ReadData
data=get(gcf,'userdata'); prm=data.prm;
[prm.td,ts]=caltomjd(prm.tstart); [tn,te]=caltomjd(prm.tend);
time=ts:prm.tint:te+(tn-prm.td)*86400; if isempty(time), return, end
prm.time=time;
data.trop={}; data.sig={}; data.ref={}; data.err={};
h=gcf; set(h,'pointer','watch');
for n=1:length(prm.rcvs)
    [data.trop{n},data.sig{n}]=...
        readtrop(prm.td,prm.time,prm.rcvs{n},prm.dirs.est,prm.fb,prm.tunit);
    if size(data.trop{n},2)==1
        data.trop{n}=[data.trop{n},zeros(length(prm.time),2)];
        data.sig{n}=[data.sig{n},repmat(nan,length(prm.time),2)];
    end
end
set(h,'name',['Tropospheric Parameters : ',prm.dirs.est,' (',prm.fb,')'])
switch prm.ref
case ''
    for n=1:length(prm.rcvs)
        data.ref{n}=repmat(nan,length(prm.time),1);
        data.err{n}=data.ref{n};
    end
otherwise
    ztd=readtrop(prm.td,prm.time,prm.rcvs,prm.dirs.trop,prm.ref);
    for n=1:length(prm.rcvs)
        data.ref{n}=ztd(:,1,n);
        data.err{n}=EstTropErr(prm.time,data.trop{n}(:,1),ztd(:,1,n));
    end
end
set(h,'pointer','arrow'); figure(h);

delete(gut('geth','rcvs'));
gut('newmenu','rcvs','&Stations',prm.rcvs,prm.rcvs,[mfilename,' cbRcv']);
gut('chkmenu','rcvs',prm.rcv);
data.prm=prm; set(gcf,'userdata',data);

% ztd estimation errors --------------------------------------------------------
function err=EstTropErr(time,est,ztd)
err=repmat(nan,length(time),1);
for i=find(~isnan(ztd))'
    j=find(time(i)-3600<=time&time<time(i)+3600); % after/before 1hrs
    if ~isempty(j), err(j)=meann(est(j))-ztd(i); end
end

% reference name ---------------------------------------------------------------
function name=RefName(ref)
names={'IGS Final','IGS Monthly','JMA MSM','JMA GSM'};
refs ={'igssnx','igsmon','mso','gso'};
i=find(strcmp(ref,refs));
if isempty(i), name=''; else name=names{i}; end

% show viewer ------------------------------------------------------------------
function Viewer(str,ti)
h=gut('newviewer','',['Summary of Tropospheric Parameters : ',ti],[600,419,0,-48],'str',str);

% plot estimated ---------------------------------------------------------------
function plotest(t,x,sig,nsig,ptype,psize)
if ptype==1|ptype==3, plot(t/3600,x,'.','markersize',psize(2)); end
if ptype==2|ptype==3, plot(t/3600,x,'-','linewidth',psize(1)); end
if nsig>0, plot(t/3600,x-nsig*sig,'r:'), plot(t/3600,x+nsig*sig,'r:'), end
range=max(abs(x))*[-1,1];
if any(isnan(range)), range=[-0.1,0.1]; end
set(gca,'ylim',range);

% plot sigma -------------------------------------------------------------------
function plotsig(t,x,sig,nsig)
if nsig>0, plot(t/3600,x-nsig*sig,'r:'), plot(t/3600,x+nsig*sig,'r:'), end

% mean/rms without nan ---------------------------------------------------------
function m=meann(x)
for n=1:size(x,2)
    i=find(~isnan(x(:,n)));
    if isempty(i), m(1,n)=nan; else m(1,n)=mean(x(i,n),1); end
end
function m=rmsn(x)
for n=1:size(x,2)
    i=find(~isnan(x(:,n)));
    if isempty(i), m(1,n)=nan; else m(1,n)=sqrt(mean(x(i,n).^2,1)); end
end
function m=mediann(x)
for n=1:size(x,2)
    i=find(~isnan(x(:,n)));
    if isempty(i), m(1,n)=nan; else m(1,n)=sqrt(median(x(i,n).^2,1)); end
end

% date/time string -------------------------------------------------------------
function str=tstr(td,ts)
str=sprintf('%04d/%02d/%02d %02d:%02d:%02.0f',mjdtocal(td,ts));

% add string -------------------------------------------------------------------
function s=adds(s,varargin), s={s{:},sprintf(varargin{:})};
