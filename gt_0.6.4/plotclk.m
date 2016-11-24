function h=plotclk(varargin)
%-------------------------------------------------------------------------------
% [system] : GpsTools
% [module] : plot satellite/receiver clocks
% [func]   : plot satellite/receiver clocks
% [argin]  : 'prm',prm   = parameters
% [argout] : h  = figure handle
% [note]   :
% [version]: $Revision: 16 $ $Date: 2008-12-12 15:49:30 +0900 (金, 12 12 2008) $
%            Copyright(c) 2004-2006 by T.Takasu, all rights reserved
% [history]: 04/11/17  0.1  new
%            04/12/15  0.2  support gui interface
%            05/06/03  0.3  add plotting allan-deviation
%            08/11/21  0.4  support igs30s,cod5s,igscod5s (gt_0.6.4)
%-------------------------------------------------------------------------------
if nargin>0&strncmp(varargin{1},'cb',2), feval(varargin{:});
else h=PlotMain(varargin{:}); end

% main -------------------------------------------------------------------------
function h=PlotMain(varargin)
h=gut('newfigg','','Satellite/Receiver Clock',[600,400,0,-68],[mfilename,' cbClose']);
if isempty(h), return, end
prm=loadprm('prm_plotclk','prm_plotclk_def');
if nargin>=2&strcmp(varargin{1},'prm'), prm=varargin{2}; end
cb1=[mfilename,' cbPlot1'];
cb2=[mfilename,' cbPlot2'];
cb3=[mfilename,' cbPlot3'];
gut('newmenu','data','&Data ',...
    {'&Read...','&New Window','-','Out&put...','&Export Plot...','-','&Options...'},{},...
    {[mfilename,' cbRead'],[mfilename,' cbNew'],'',[mfilename,' cbOut'],...
     [mfilename,' cbExport'],'',[mfilename,' cbOpt']});
gut('newmenu','plot','&Plot',...
    {'Clock &Bias','Clock &Drift','Clock Bias &Error','-','Clock Error by &Sats/Rcvs',...
     'Clock Error by &Dates','Clock S&tability','-','Su&mmary...'},...
    {'clkbias','clkdrift','clkerr','-','errsat','errdate','allandev','',''},...
    {cb1,cb1,cb1,'',cb2,cb2,cb2,'',cb3});
gut('newmenu','sats','&Satellite',prm.sats,prm.sats,[mfilename,' cbSat']);
gut('newmenu','rcvs','&Receiver',prm.rcvs,prm.rcvs,[mfilename,' cbRcv']);
gut('chkmenu','plot',prm.type);
gut('chkmenu','sats',prm.sat);
gut('chkmenu','rcvs',prm.rcv);
if strcmp(prm.sr,'sat'), gut('hidemenu','rcvs'); else gut('hidemenu','sats'); end
set(gut('newbtnh','sbtn',[0,0,36,16],{'<','>'},[mfilename,' cbSwBtn']),...
    'fontsize',8);
set(gcf,'resizefcn',[mfilename,' cbResize']);
LocatePos;
data.clk={}; data.err={};
data.prm=prm; set(h,'userdata',data);
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
function cbResize, LocatePos

% callback on menu read --------------------------------------------------------
function cbRead
sel1={{'Satellite Clock','Receiver Clock'},{'sat','rcv'}};
sel2={{'OFF','ON'},0:1};
sel3={{'','Ref Rcv-Clock','Mean Sat-Clock'},{'','rcv','sat'}};
e1={'Estimated(Forward)','Estimated(Backward)','Estimated(Smoothed)'};
e2={'clkf','clkb','clkfb'};
r1={'IGS Final','IGS Rapid','IGS URapid','IGS URapid(Pred)','CODE','EMR','ESA',...
    'GFZ','JPL','MIT','NGS','IGS/CODE','CODE Rapid','IGR/CODR','IGS 30s',...
    'CODE 5s','IGS/CODE 5s','Broadcast'};
r2={'igs','igr','igu','igp','cod','emr','esa',...
    'gfz','jpl','mit','ngs','igscod','codr','igrcodr','igs30s',...
    'cod5s','igscod5s','brdc'};
data=get(gcf,'userdata'); prm=data.prm;
prm1={
't','Start Time (GPST)',[]
't','End Time (GPST)',[]
'n','Time Interval (sec)', 0
'n','Processing Unit Time (hr)', 0
'p','Clock Type', sel1
'b','Satellites/Receivers',[mfilename,' cbSelSat']
'p','Satellite/Receiver Clock',{{e1{:},r1{:}},{e2{:},r2{:}}}
'p','Reference Satellite/Receiver Clock',{{'',e1{:},r1{:}},{'',e2{:},r2{:}}}
'p','Baseline Clock', sel3
'p','Interpolate Reference Clock', sel2
};
prm2={
' ','Satellite/Receiver Clock Directory',''
'd','',prm.dirs.est
' ','Reference Satellite/Receiver Clock Directory',''
'd','',prm.dirs.clk
};
q='';
gut('newdlg','','Read Data',[320,366]);
gut('newprms','prm1',[12,112,300,25,120],prm1);
gut('newprms','prm2',[12,33,300,18,120],prm2);
gut('newokcancelbtn','',[132,4,180,23]);
gut('setprms','prm1',prm.tstart,prm.tend,prm.tint,prm.tunit,prm.sr,'',prm.fb,prm.ref,prm.refclk,prm.interp);
gut('setudata','prm1_6',{prm.sats,prm.rcvs});
if ~gut('waitok'), return, end
[prm.tstart,prm.tend,prm.tint,prm.tunit,prm.sr,q,prm.fb,prm.ref,prm.refclk,prm.interp]=gut('getprms','prm1');
[q,prm.dirs.est,q,prm.dirs.clk]=gut('getprms','prm2');
list=gut('getudata','prm1_6'); prm.sats=list{1}; prm.rcvs=list{2};
if ~isempty(prm.sats), prm.sat=prm.sats{1}; else prm.sat=''; end
if ~isempty(prm.rcvs), prm.rcv=prm.rcvs{1}; else prm.rcv=''; end
close
switch prm.sr
case 'sat', gut('hidemenu','rcvs'); gut('showmenu','sats');
case 'rcv', gut('hidemenu','sats'); gut('showmenu','rcvs');
end
data.prm=prm; set(gcf,'userdata',data);
SaveSetting; ReadData; UpdatePlot;

% callback on menu select satellite ------------------------------------------
function cbSelSat
list=get(gcbo,'userdata');
switch gut('getsel','prm1_5')
case 'sat', [list{1},ok]=editlist('Satellites',list{1});
case 'rcv', [list{2},ok]=editlist('Receivers',list{2},'rcv');
end
if ok, set(gcbo,'userdata',list); end

% callback on menu output results ----------------------------------------------
function cbOut
C=299792458;
data=get(gcf','userdata'); prm=data.prm;
if strcmp(prm.sr,'sat')
    i=find(strcmp(prm.sat,prm.sats)); ti='SATELLITE'; name=prm.sats{i};
else
    i=find(strcmp(prm.rcv,prm.rcvs)); ti='RECEIVER'; name=prm.rcvs{i};
end
if isempty(data.clk)|isempty(i), return, end
s1={};
s1=adds(s1,'%%  %s CLOCK : %s : %s-%s GPST',ti,name,tstr(prm.td,prm.time(1)),tstr(prm.td,prm.time(end)));
s1=adds(s1,'%%  date     time(sec)   clock-bias (sec)       std (sec)');
s1=adds(s1,'%%---------------------------------------------------------');
s2=cell(1,length(prm.time));
for n=1:length(prm.time)
    d=mjdtocal(prm.td,prm.time(n)); t=mod(prm.time(n),86400);
    s2{n}=sprintf('%04d,%02d,%02d,%9.2f, %17.13f, %17.13f\n',d(1:3),t,...
                  data.clk{i}(n,1)/C,data.sig{i}(n,1)/C);
end
Viewer({s1{:},s2{:}},['Satellite/Receiver Clock : ',prm.dirs.est,' (',prm.fb,')']);

% callback on menu export plot -------------------------------------------------
function cbExport
persistent file dirs
if isempty(file), file='plotclk'; end, if isempty(dirs), dirs=''; end
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
sr=strcmp(prm.sr,'sat');
if isempty(strmatch(prm.type,'clk'))
    if sr, sat=prm.sat; m=length(prm.sats); else rcv=prm.rcv; m=length(prm.rcvs); end
    for n=1:m
        if sr, prm.sat=prm.sats{n}; else prm.rcv=prm.rcvs{n}; end
        data.prm=prm; set(f,'userdata',data)
        figure(f), UpdatePlot;
        print(format{:},opts{:},sprintf('%s%02d',fs,n));
    end
    if sr, prm.sat=sat; else prm.rcv=rcv; end
    data.prm=prm; set(f,'userdata',data); UpdatePlot;
else
    print(format{:},opts{:},fs);
end

% callback on menu options -----------------------------------------------------
function cbOpt
C=299792458;
sel1={{'OFF','ON'},0:1};
sel3={{'m','nsec'},0:1};
sel4={{'Dot','Line','Line and Dot'},1:3};
data=get(gcf,'userdata'); prm=data.prm;
prm1={
'p','Show Reference Clock',            sel1
'p','Show Mean and RMS Error',         sel1
'p','Show Relative Values',            sel1
'n','Show Est. Std Dev. (sigma,0:OFF)',0
'p','Clock Units',                     sel3
'p','Plot Style',                      sel4
's','Plot Line Width / Marker Size',   [0.1,10]
's','Plot Axis Range (nsec)',          [-1,1]
};
gut('newdlg','','Options',[320,238]);
gut('newprms','prm1',[12,34,300,25,120],prm1);
gut('newokcancelbtn','',[133,4,180,23]);
gut('setprms','prm1',prm.showf(1),prm.showf(2),prm.showf(3),prm.nsig,prm.sec,...
    prm.ptype,prm.psize,prm.range)
if ~gut('waitok'), return, end
[prm.showf(1),prm.showf(2),prm.showf(3),prm.nsig,prm.sec,prm.ptype,prm.psize,...
 range]=gut('getprms','prm1');
if range(1)<range(2), prm.range=range; else prm.range=[-1,1]; end
close
data.prm=prm; set(gcf,'userdata',data);
SaveSetting; UpdatePlot;

% callback on menu plot --------------------------------------------------------
function cbPlot1
types={'errsat','errdate','allandev'};
data=get(gcf,'userdata'); prm=data.prm;
gut('unchkmenu','plot',types);
gut('togglechk',gcbo);
prm.type=gut('getchk','plot');
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

% callback on menu satellite ---------------------------------------------------
function cbSat
data=get(gcf,'userdata'); prm=data.prm;
prm.sat=get(gcbo,'userdata');
gut('unchkmenu','sats');
gut('chkmenu','sats',prm.sat);
data.prm=prm; set(gcf,'userdata',data); UpdatePlot;

% callback on menu receiver ----------------------------------------------------
function cbRcv
data=get(gcf,'userdata'); prm=data.prm;
prm.rcv=get(gcbo,'userdata');
gut('unchkmenu','rcvs');
gut('chkmenu','rcvs',prm.rcv);
data.prm=prm; set(gcf,'userdata',data); UpdatePlot;

function EditRcv
data=get(gcf,'userdata'); prm=data.prm;
[prm.rcvs,ok]=editlist('Receivers',prm.rcvs,'rcv');
if ~ok, return, end
delete(gut('geth','rcvs'));
gut('newmenu','rcvs','&Receivers',...
    {prm.rcvs{:},'-','&Edit...'},{prm.rcvs{:},'',''},[mfilename,' cbRcv']);
gut('chkmenu','rcvs',prm.rcv);
data.prm=prm; set(gcf,'userdata',data); ReadData; UpdatePlot;

% callback on push swithch button ----------------------------------------------
function cbSwBtn
data=get(gcf,'userdata'); prm=data.prm;
if strcmp(prm.sr,'sat')
    j=find(strcmp(prm.sats,prm.sat));
    switch get(gcbo,'tag')
    case 'sbtn_1', prm.sat=prm.sats{max(j-1,1)};
    case 'sbtn_2', prm.sat=prm.sats{min(j+1,length(prm.sats))};
    end
    gut('unchkmenu','rcvs'); gut('chkmenu','rcvs',prm.rcv);
else
    j=find(strcmp(prm.rcvs,prm.rcv));
    switch get(gcbo,'tag')
    case 'sbtn_1', prm.rcv=prm.rcvs{max(j-1,1)};
    case 'sbtn_2', prm.rcv=prm.rcvs{min(j+1,length(prm.rcvs))};
    end
    gut('unchkmenu','rcvs'); gut('chkmenu','rcvs',prm.rcv);
end
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
save(fullfile(path,'settings','prm_plotclk.mat'),'prm');

% update plot ------------------------------------------------------------------
function UpdatePlot
data=get(gcf,'userdata'); prm=data.prm;
if any(strcmp(prm.type,'errsat')), PlotEstErrSat
elseif any(strcmp(prm.type,'errdate')), PlotEstErrDate
elseif any(strcmp(prm.type,'allandev')), PlotAllanDev
else PlotClock; end

% plot clock -------------------------------------------------------------------
function PlotClock
data=get(gcf,'userdata'); prm=data.prm;
if isempty(data.clk), return, end
clk=[]; sig=[]; err=[]; ref=[];
ti=[tstr(prm.td,prm.time(1)),'-',tstr(prm.td,prm.time(end))];
if strcmp(prm.sr,'sat')
     ti=['Satellite Clock ',prm.sat,' : ',ti]; i=find(strcmp(prm.sat,prm.sats));
else ti=['Receiver Clock ',prm.rcv,' : ',ti]; i=find(strcmp(prm.rcv,prm.rcvs)); end
if ~isempty(i)
    clk=data.clk{i}; sig=data.sig{i}; err=data.err{i}; rclk=data.rclk{i};
end
clf
if isempty(clk), return, end
C=299792458;
margin=[0.08,0.03,0.04,0.012];
if prm.sec, clk=clk*1E9/C; sig=sig*1E9/C; err=err*1E9/C; rclk=rclk*1E9/C; end
for n=1:length(prm.type)
    if n<length(prm.type), topts='nolabel'; else topts=''; end
    ggt('subplotv','est',n,length(prm.type),'move',[1,1,1,1],'link',[1,0,1,0],...
        'taxis',prm.td,'topts',topts,'margin',margin,...
        'xlim',[prm.time(1),prm.time(end)]/3600);
    switch prm.type{n}
    case 'clkbias'
        if prm.showf(3), off=meann(clk(:,1)); else off=0; end
        range=[min(clk(:,1)),max(clk(:,1))]-off;
        if prm.showf(1)&~isempty(prm.ref)
            h=plot(prm.time/3600,rclk-off,'m.','markersize',prm.psize(2));
            legend(h,RefName(prm.ref,prm.interp));
            range=[min([range(1);rclk-off]),max([range(2);rclk-off])];
        end
        if prm.ptype==1|prm.ptype==3, plot(prm.time/3600,clk(:,1)-off,'.','markersize',prm.psize(2)); end
        if prm.ptype==2|prm.ptype==3, plot(prm.time/3600,clk(:,1)-off,'-','linewidth',prm.psize(1)); end
        if any(isnan(range)), range=[0,1]; end
        set(gca,'ylim',range);
        if prm.sec, unit='nsec'; else unit='m'; end
        yl='Clock Bias'; if prm.showf(3), yl=sprintf('%s-%.4f',yl,off); end
        ylabel([yl,' (',unit,')'])
    case 'clkdrift'
        if prm.showf(3), off=meann(clk(:,2)); else off=0; end
        if size(clk,2)>=2
            if prm.ptype==1|prm.ptype==3, plot(prm.time/3600,clk(:,2)-off,'.','markersize',prm.psize(2)); end
            if prm.ptype==2|prm.ptype==3, plot(prm.time/3600,clk(:,2)-off,'-','linewidth',prm.psize(1)); end
            range=[min(clk(:,2)-off),max(clk(:,2)-off)];
            if any(isnan(range)), range=[0,1]; end
            set(gca,'ylim',range);
        end
        if prm.sec, unit='nsec/sec'; else unit='m/sec'; end
        yl='Clock Drifts'; if prm.showf(3), yl=sprintf('%s-%.4f',yl,off); end
        ylabel([yl,' (',unit,')'])
    case 'clkerr'
        if prm.showf(3), off=meann(err(:,1)); else off=0; end
        if prm.ptype==1|prm.ptype==3, plot(prm.time/3600,err(:,1)-off,'.','markersize',prm.psize(2)); end
        if prm.ptype==2|prm.ptype==3, plot(prm.time/3600,err(:,1)-off,'-','linewidth',prm.psize(1)); end
        plotsig(prm.time,0,sig(:,1),prm.nsig)
        if prm.sec, range=prm.range; unit='nsec'; else range=prm.range/1E9*C; unit='m'; end
        if isinf(range(1)), range(1)=max(err(:,1)-off); end
        if isinf(range(2)), range(2)=min(err(:,1)-off); end
        set(gca,'ylim',range);
        yl='Clock Bias Error'; if prm.showf(3), yl=sprintf('%s-%.4f',yl,off); end
        ylabel([yl,' (',unit,')'])
        if prm.showf(2)
            s=sprintf('REF: %s  MEAN: %6.4f%s RMS: %6.4f%s',RefName(prm.ref,prm.interp),...
                      meann(err),unit,rmsn(err),unit);
            ggt('mtext',['est_',num2str(n)],s,1);
        end
    end
    if n==1, title(ti); end
end

% plot estimation error by satellites ------------------------------------------
function PlotEstErrSat
data=get(gcf,'userdata'); prm=data.prm;
C=299792458; units='(nsec)';
if strcmp(prm.sr,'sat'), sr='Satellite'; names=prm.sats; 
else sr='Receiver'; names=prm.rcvs; end 
pos=[0.08,0.12,0.89,0.81];
clf
if isempty(data.err), return, end
for n=1:length(data.err) err(n,:)=rmsn(data.err{n}(:,1)); end
err=[err;meann(err)]; range=prm.range-prm.range(1); 
if prm.sec, err=err/C*1E9; else range=range*C/1E9; units='(m)'; end
h=ggt('barplot',err,{names{:},'Average'},'ylim',range,'position',pos);
legend(h,'Clock Bias')
ylabel(['Clock Bias RMS Error ',units]);
ti=[sr,' Clock : ',tstr(prm.td,prm.time(1)),'-',tstr(prm.td,prm.time(end)),...
    ' (REF: ',RefName(prm.ref,prm.interp),')'];
title(ti);

% plot allan deviation ---------------------------------------------------------
function PlotAllanDev
data=get(gcf,'userdata'); prm=data.prm;
C=299792458;
if strcmp(prm.sr,'sat'), sr='Satellite'; names=prm.sats; 
else sr='Receiver'; names=prm.rcvs; end 
tau=[1,3,10,30,100,300,600,1000,3000,6000,10000,30000,60000,100000,300000];
sig=repmat(nan,length(tau),length(names));
for n=1:length(names)
    for m=1:length(tau)
        if tau(m)>=prm.tint&mod(tau(m),prm.tint)==0
            i=1:tau(m)/prm.tint:length(prm.time);
            if length(i)>1
                clk=data.clk{n}(i,1)/C;
                sig(m,n)=stdn((clk(1:end-1)-clk(2:end))/tau(m))/sqrt(2);
            end
        end
    end
end
clf
t=repmat(tau',1,length(names));
loglog(t,sig,'-','linewidth',prm.psize(1)); hold on
h=loglog(t,sig,'.','markersize',prm.psize(2));
legend(h,names);
[fn,fs]=gut('getfont','g');
set(gca,'position',[0.09,0.1,0.88,0.85],'fontname',fn,'fontsize',fs);
xlim([1E-1,1E5]), ylim([1E-15,1E-10]);
grid on
xlabel('Averaging Time \tau (sec)');
ylabel('Allan Deviation \sigma_y (\tau)');
ti=[sr,' Clock : ',tstr(prm.td,prm.time(1)),'-',tstr(prm.td,prm.time(end))];
title(['Stability of ',ti]);

% plot estimatin error by date -------------------------------------------------
function PlotEstErrDate
data=get(gcf,'userdata'); prm=data.prm;
C=299792458; units='(nsec)';
if strcmp(prm.sr,'sat'), sr='Satellite'; else sr='Receiver'; end 
pos=[0.08,0.12,0.89,0.81];
clf
if isempty(data.err), return, end
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
err=[err;meann(err)]; range=prm.range-prm.range(1); 
if prm.sec, err=err/C*1E9; else range=range*C/1E9; units='(m)'; end
h=ggt('barplot',err,{labels{:},'Average'},'ylim',range,'position',pos);
legend(h,'Clock Bias')
ylabel(['Clock Bias RMS Error ',units]);
ti=[sr,' Clock : ',tstr(prm.td,prm.time(1)),'-',tstr(prm.td,prm.time(end)),...
    ' (REF: ',RefName(prm.ref,prm.interp),')'];
title(ti);

% update summary ---------------------------------------------------------------
function UpdateSummary
data=get(gcf','userdata'); prm=data.prm;
C=299792458;
if isempty(data.err), return, end
if strcmp(prm.sr,'sat'), type='SATELLITE'; else type='RECEIVER'; end
s={};
s=adds(s,'                %s CLOCK : CLOCK BIAS ERROR (REF: %s)',type,RefName(prm.ref,prm.interp));
s=adds(s,'');
s=AddSummary(s,data,prm,prm.time(1),prm.time(end),'');
s=adds(s,'');
s=adds(s,' -----------------------------------------------------------------------------');
tu=prm.tunit*3600;
ts=(floor(prm.time(1)/tu):floor(prm.time(end)/tu))*tu;
for n=1:length(ts)
    s=adds(s,'');
    if isempty(data.estref), estref='GPST'; else estref=data.estref{n}; end
    s=AddSummary(s,data,prm,ts(n),ts(n)+tu-prm.tint,estref);
end
ti=[prm.dirs.est,' (',prm.fb,')'];
Viewer(s,['Summary of Satellite/Receiver Clock : ',ti]);

% add summary string -----------------------------------------------------------
function s=AddSummary(s,data,prm,ts,te,estref)
C=299792458;
if strcmp(prm.sr,'sat'), type='SATELLITE'; names=prm.sats;
else type='RECEIVER'; names=prm.rcvs; end
if prm.sec, unit='(ns)'; f=1E9/C; else unit='(m) '; f=1; end
s=adds(s,' TIME    : %s-%s GPST',tstr(prm.td,ts),tstr(prm.td,te));
if ~isempty(estref), s=adds(s,' REFCLK   : %s',estref); end
s=adds(s,'              mean%s  rms%s  mean*%s  rms*%s  est.rate   ref.rate',unit,unit,unit,unit);
s=adds(s,' -----------------------------------------------------------------------------');
p=[]; v=[]; pa=[]; pv=[]; r=[];
i=find(ts<=prm.time&prm.time<=te);
for n=1:length(names)
   p(n,1)=meann(data.err{n}(i,1))*f;
   v(n,1)=rmsn(data.err{n}(i,1))*f;
   r(n,1)=length(find(~isnan(data.clk{n}(i,1))))/length(i)*100;
   rr(n,1)=length(find(~isnan(data.rclk{n}(i,1))))/length(i)*100;
end
for n=1:length(names)
   mp=meann(data.err{n}(i,1));
   pa(n,1)=meann(data.err{n}(i,1)-mp)*f;
   va(n,1)=rmsn(data.err{n}(i,1)-mp)*f;
   if strcmp(names{n},data.estref)
       c='refclk';
       s=adds(s,' %-7s : %9s %9s %9s %9s %10s %9.1f%%',names{n},c,c,c,c,c,rr(n));
   else
       s=adds(s,' %-7s : %9.4f %9.4f %9.4f %9.4f %9.1f%% %9.1f%%',names{n},...
              p(n),v(n),pa(n),va(n),r(n),rr(n));
   end
end
s=adds(s,' -----------------------------------------------------------------------------');
s=adds(s,' %-7s : %9.4f %9.4f %9.4f %9.4f %9.1f%% %9.1f%%','average',...
       meann(p),meann(v),meann(pa),meann(va),meann(r),meann(rr));
s=adds(s,' %-7s : %9.4f %9.4f %9.4f %9.4f %9.1f%% %9.1f%%','median',...
       mediann(p),mediann(v),mediann(pa),mediann(va),mediann(r),meann(rr));

% read clock data --------------------------------------------------------------
function ReadData
data=get(gcf,'userdata'); prm=data.prm;
[prm.td,ts]=caltomjd(prm.tstart); [tn,te]=caltomjd(prm.tend);
time=ts:prm.tint:te+(tn-prm.td)*86400; if isempty(time), return, end
prm.time=time;
gut('newmsgbar','msgbar','',[480,50],[0,0.7,1]);
gut('setmsgbar','msgbar','reading satellite/receiver clock',0);
if strcmp(prm.sr,'sat'), names=prm.sats; else names=prm.rcvs; end
[data.clk,data.sig,data.err,data.rclk,data.estref]=...
    ReadEst(prm.td,prm.time,names,prm.dirs,prm.fb,prm.tunit,prm.sr,prm.ref,...
            prm.refclk,prm.interp);

gut('closemsgbar','msgbar');
data.prm=prm; set(gcf,'userdata',data);

% update menu
if strcmp(prm.sr,'sat')
    delete(gut('geth','sats'));
    gut('newmenu','sats','&Satellite',prm.sats,prm.sats,[mfilename,' cbSat']);
    gut('chkmenu','sats',prm.sat);
    ti=['Satellite Clock : ',prm.dirs.est,' (',prm.fb,')'];
else
    delete(gut('geth','rcvs'));
    gut('newmenu','rcvs','&Receiver',prm.rcvs,prm.rcvs,[mfilename,' cbRcv']);
    gut('chkmenu','rcvs',prm.rcv);
    ti=['Receiver Clock : ',prm.dirs.est,' (',prm.fb,')'];
end
set(gcf,'name',ti);

% read estimated clock ---------------------------------------------------------
function [clk,sig,err,rclk,refrcv]=ReadEst(td,time,names,dirs,fb,tunit,sr,ref,...
                                           refclk,interp)
clk={}; sig={}; err={}; rclk={}; refrcv={}; tu=tunit*3600;
ts=(floor(time(1)/tu):floor(time(end)/tu))*tu;
for n=1:length(ts)
    i=find(ts(n)<=time&time<ts(n)+tu);
    [clkn,sign,rclkn,errn,refn]=...
        ReadEstUnit(td,time(i),names,dirs,fb,tunit,sr,ref,refclk,interp);
    for m=1:length(names)
        clk{m}(i,:)=clkn(:,:,m);
        sig{m}(i,:)=sign(:,:,m);
        err{m}(i,:)=errn(:,:,m);
        rclk{m}(i,:)=rclkn(:,:,m);
    end
    refrcv={refrcv{:},refn};
end

% read clock data in processing unit time --------------------------------------
function [clk,sig,rclk,err,refrcv]=ReadEstUnit(td,time,names,dirs,fb,tunit,sr,...
                                               ref,refclk,interp)
C=299792458;
[clk,sig,refrcv]=readclk(td,time,names,dirs.est,fb,tunit);
clk=C*clk; sig=C*sig;

if interp, opt='interp'; else opt=''; end

if isempty(ref)
    rclk=repmat(nan,size(clk)); err=rclk;

elseif strcmp(refclk,'rcv')&~isempty(refrcv) % relatve to reference clock station
    rclk=C*readclk(td,time,{names{:},refrcv},dirs.clk,ref,tunit,opt);
    for n=1:length(names), clk(:,1,n)=clk(:,1,n)+rclk(:,:,end); end
    rclk=rclk(:,:,1:end-1);
    err=clk(:,1,:)-rclk(:,1,:);
elseif strcmp(refclk,'sat')&strcmp(sr,'sat') % relative to satclock average
    rclk=C*readclk(td,time,names,dirs.clk,ref,tunit,opt);
    err=clk-rclk;
    for n=1:size(err,1), err(n,:,:)=err(n,:,:)-meann(squeeze(err(n,:,:))); end
else
    rclk=C*readclk(td,time,names,dirs.clk,ref,tunit,opt);
    err=clk(:,1,:)-rclk(:,1,:);
end

% reference name ---------------------------------------------------------------
function name=RefName(ref,interp)
names={'IGS Final','IGS Rapid','IGS URapid','JPL','IGS-COD'};
refs ={'igs','igr','igu','jpl','igscod'};
i=find(strcmp(ref,refs));
if isempty(i), name=ref; else name=names{i}; end
if interp, name=[name,'(interpolated)']; end

% show viewer ------------------------------------------------------------------
function Viewer(str,ti)
h=gut('newviewer','',ti,[600,419,0,-48],'str',str);

% plot sigma -------------------------------------------------------------------
function plotsig(t,x,sig,nsig)
if nsig>0, plot(t/3600,x-nsig*sig,'r:'), plot(t/3600,x+nsig*sig,'r:'), end

% mean/rms without nan ---------------------------------------------------------
function m=meann(x)
for n=1:size(x,2)
    i=find(~isnan(x(:,n)));
    if isempty(i), m(1,n)=nan; else m(1,n)=mean(x(i,n),1); end
end
function m=stdn(x)
for n=1:size(x,2)
    i=find(~isnan(x(:,n)));
    if isempty(i), m(1,n)=nan; else m(1,n)=std(x(i,n),1); end
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
