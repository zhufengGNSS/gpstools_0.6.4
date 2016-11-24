function h=plotpos(varargin)
%-------------------------------------------------------------------------------
% [system] : GpsTools
% [module] : plot receiver positions
% [func]   : plot receiver positionns
% [argin]  : 'prm',prm   = parameters
% [argout] : h = figure handle
% [note]   :
% [version]: $Revision: 31 $ $Date: 06/07/20 19:10 $
%            Copyright(c) 2004-2006 by T.Takasu, all rights reserved
% [history]: 04/11/17  0.1  new
%-------------------------------------------------------------------------------
if nargin>0&strncmp(varargin{1},'cb',2), feval(varargin{:})
else h=PlotMain(varargin{:}); end

% main -------------------------------------------------------------------------
function h=PlotMain(varargin)
h=gut('newfigg','','Receiver Position',[600,400,0,-68],[mfilename,' cbClose']);
if isempty(h), return, end
prm=loadprm('prm_plotpos','prm_plotpos_def');
if nargin>=2&strcmp(varargin{1},'prm'), prm=varargin{2}; end
cb1=[mfilename,' cbPlot1'];
cb2=[mfilename,' cbPlot2'];
cb3=[mfilename,' cbPlot3'];
gut('newmenu','data','&Data ',...
    {'&Read...','&New Window','-','Out&put...','&Export Plot...','-',...
     '&Map Area...','&Options...'},{},...
    {[mfilename,' cbRead'],[mfilename,' cbNew'],'',[mfilename,' cbOut'],...
     [mfilename,' cbExport'],'',[mfilename,' cbMap'],[mfilename,' cbOpt']});
gut('newmenu','plot','&Plot',...
    {'&Position Error','Position Error &X/Y/Z','Position Error &E/N/U',...
     'Position &Displacement','Position Displacement &3D',...
     'Position Error &East','Position Error &North','Position Error &Up','-',...
     'Position Mean Error by Receivers','Position RMS Error by &Receivers',...
     'Position Error by &Dates','-','&Horizontal Error Map','&Vertical Error Map',...
     'Horizontal &Variation Map','Vertical V&ariation Map','&Trans. Parameters',...
     '-','&Summary...','Transformation Params...'},...
    {'posh','posx','post','posw','pos3','pose','posn','posu','','posm','posr',...
     'posd','','mapeh','mapev','maph','mapv','trpp','','',''},...
    {cb1,cb1,cb1,cb1,cb1,cb1,cb1,cb1,'',cb1,cb1,cb1,'',cb1,cb1,cb1,cb1,cb1,'',...
     cb2,cb3});
gut('newmenu','rcvs','&Receivers',prm.rcvs,prm.rcvs,[mfilename,' cbRcv']);
gut('chkmenu','plot',prm.type);
gut('chkmenu','rcvs',prm.rcv);
set(gut('newbtnh','sbtn',[0,0,36,16],{'<','>'},[mfilename,' cbSwBtn']),...
    'fontsize',8);
set(gcf,'resizefcn',[mfilename,' cbResize']);
LocatePos;
data.pos={}; data.err={}; data.sig={}; data.posr={}; data.timef=[]; data.tr=[];
data.prm=prm; set(gcf,'userdata',data);
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
sel1={{'Estimated(Forward)','Estimated(Backward)','Estimated(Smoothed)',...
       'EstSinex(Forward)','EstSinex(Backward)','EstSinex(Smoothed)',...
       'ITRF97','ITRF2000','IGS00','IGb00','IGS Final','GSIPOS'},...
      {'posf','posb','posfb','possf','possb','possfb',...
       'itrf97','itrf2000','igs00','igb00','igssnx','gsipos'}};
sel2={{'Average','Linear Fit','Poly Fit','Start','End','Est Average','ITRF97','ITRF2000',...
       'IGS00','IGb00','IGS Final','GSIPOS'},...
      {'estm','estf','estp','ests','este','est0','itrf97','itrf2000','igs00','igb00','igssnx','gsipos'}};
sel3={{'','ITRF97','IGS00','IGb00'},{'','itrf97','igs00','igb00'}};
sel4={{'OFF','ON'},0:1};
data=get(gcf,'userdata'); prm=data.prm;
prm1={
't','Start Time (GPST)',[]
't','End Time (GPST)',[]
'n','Time Interval (sec)',0
'n','Processing Unit Time (hr)',0
'b','Receivers',[mfilename,' cbSelRcv']
'p','Receiver/Station Postions',sel1
'p','Reference Postions',sel2
'p','Adjust To Reference Frame',sel3
's','Bandpass Filter (Hz)',[0,0]
's','Sidereal Filter Days/Ajust (day,sec)',[0,0,0]
's','Sidereal Filter Bandpass (Hz)', [0,0]
};
prm2={
' ','Receiver Position Directory',''
'd','',prm.dirs.est
' ','Reference Station Position Directory',''
'd','',prm.dirs.pos
' ','Coordinate Transformation Parameters',''
'h','',prm.trprm
};
gut('newdlg','','Read Data',[314,400]);
gut('newprms','prm1',[12,140,295,23,118],prm1);
gut('newprms','prm2',[12,32,295,17],prm2);
gut('newokcancelbtn','',[127,4,180,23]);
gut('setprms','prm1',prm.tstart,prm.tend,prm.tint,prm.tunit,'',prm.fb,prm.ref,prm.adj,prm.bpfilt,prm.srday,prm.srprm);
gut('setudata','prm1_5',prm.rcvs);
if ~gut('waitok'), return, end
[prm.tstart,prm.tend,prm.tint,prm.tunit,q,prm.fb,prm.ref,prm.adj,prm.bpfilt,prm.srday,prm.srprm]=gut('getprms','prm1');
[q,prm.dirs.est,q,prm.dirs.pos,q,prm.trprm]=gut('getprms','prm2');
prm.rcvs=gut('getudata','prm1_5');
if isempty(prm.rcvs), prm.rcv=''; else prm.rcv=prm.rcvs{1}; end
prm.files={};
close
data.prm=prm; set(gcf,'userdata',data);
SaveSetting; ReadData; UpdatePlot;

% callback on menu select station ---------------------------------------------
function cbSelRcv
[rcvs,ok]=editlist('Receivers',get(gcbo,'userdata'),'rcv');
if ok, set(gcbo,'userdata',rcvs); end

% output position --------------------------------------------------------------
function cbOut
data=get(gcf','userdata'); prm=data.prm;
i=find(strcmp(prm.rcv,prm.rcvs));
if isempty(data.pos)|isempty(i), return, end
s1={};
s1=adds(s1,'%%  RECEIVER POSITION : %s : %s-%s GPST',RcvName(prm.rcvs{i},prm.showf(4)),tstr(prm.td,prm.time(1)),tstr(prm.td,prm.time(end)));
s1=adds(s1,'%%  date     time(sec)       x(m)           y(m)           z(m)       sdx(m)  sdy(m)  sdz(m)');
s1=adds(s1,'%%------------------------------------------------------------------------------------------');
s2=cell(1,length(prm.time));
for n=1:length(prm.time)
    d=mjdtocal(prm.td,prm.time(n)); t=mod(prm.time(n),86400);
    s2{n}=sprintf('%04d,%02d,%02d,%9.2f, %14.4f,%14.4f,%14.4f, %7.4f,%7.4f,%7.4f',...
              d(1:3),t,data.pos{i}(n,1:3),data.sig{i}(n,1:3));
end
Viewer({s1{:},s2{:}},['Receiver Position : ',prm.dirs.est,' (',prm.fb,')']);

% callback on menu export plot -------------------------------------------------
function cbExport
persistent file dirs
if isempty(file), file='plotpos'; end, if isempty(dirs), dirs=''; end
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
if any(strcmp(prm.type,{'posh','posx','post','posw','pos3'}))
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
data=get(gcf,'userdata'); prm=data.prm;
sel1={{'OFF','ON'},0:1};
sel2={{'Dot','Line','Line and Dot','Circle'},1:4};
sel4={{'OFF','Code Name','Full Name','Both Name'},0:3};
sel5={{'All','Final Only'},0:1};
prm1={
'p','Show Estimated Position',sel5
'p','Show Reference Position',sel1
'p','Show Mean and RMS Error',sel1
'p','Show Estimation Date',sel1
'p','Show Receiver Name in Map',sel4
'n','Show Est. Std Dev. (sigma,0:OFF)',0
'p','Plot Style',sel2
's','Plot Line Width / Marker Size',[0.5,5]
'g','Receiver Name Font in Map',''
'p','Detrend Position',sel1
's','Plot Axis Range',[-0.04,0.04]
};
gut('newdlg','','Options',[320,295]);
gut('newprms','prm1',[12,35,300,23,120],prm1);
gut('newokcancelbtn','',[133,4,180,23]);
gut('setprms','prm1',prm.est,prm.showf(1),prm.showf(2),prm.showf(3),prm.showf(4),...
    prm.nsig,prm.ptype,prm.psize,prm.rfont,prm.detr,prm.range)
if ~gut('waitok'), return, end
[prm.est,prm.showf(1),prm.showf(2),prm.showf(3),prm.showf(4),prm.nsig,prm.ptype,...
 prm.psize,prm.rfont,prm.detr,prm.range]=gut('getprms','prm1');
close
data.prm=prm; set(gcf,'userdata',data);
SaveSetting; UpdatePlot;

% callback on menu map area ----------------------------------------------------
function cbMap
data=get(gcf,'userdata'); prm=data.prm;
[prm.map,ok]=editmap('Map Area',prm.map);
if ok, data.prm=prm; set(gcf,'userdata',data); UpdatePlot; end

% callback on menu plot --------------------------------------------------------
function cbPlot1
data=get(gcf,'userdata'); prm=data.prm;
gut('unchkmenu','plot');
set(gcbo,'checked','on');
prm.type=get(gcbo,'userdata');
data.prm=prm; set(gcf,'userdata',data); UpdatePlot;

function cbPlot2, SummaryPosErr
function cbPlot3, OutputTrp

% callback on menu station -----------------------------------------------------
function cbRcv
data=get(gcf,'userdata'); prm=data.prm;
prm.rcv=get(gcbo,'userdata');
if isempty(prm.rcv), EditRcv, return, end
gut('unchkmenu','rcvs'); gut('chkmenu','rcvs',prm.rcv);
data.prm=prm; set(gcf,'userdata',data); UpdatePlot;

function EditRcv
data=get(gcf,'userdata'); prm=data.prm;
[prm.rcvs,ok]=editlist('Receivers',prm.rcvs,'rcv');
if ok, data.prm=prm; set(gcf,'userdata',data); UpdatePlot; end

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
UpdatePlot;

% locate position --------------------------------------------------------------
function LocatePos
p=get(gcf,'position'); m=15;
gut('setpos','sbtn_1',[p(3)-m-32,p(4)-16]);
gut('setpos','sbtn_2',[p(3)-m-16,p(4)-16]);

% save settingsrs --------------------------------------------------------------
function SaveSetting
data=get(gcf,'userdata'); prm=data.prm;
[path,file]=fileparts(which(mfilename));
save(fullfile(path,'settings','prm_plotpos.mat'),'prm');

% update plot/display ----------------------------------------------------------
function UpdatePlot
data=get(gcf','userdata');
switch data.prm.type
case 'posh', PlotPosH
case 'pos3', PlotPos3D
case 'posx', PlotPosX
case 'post', PlotPosE
case 'posw', PlotPosW
case 'pose', PlotPosS(1)
case 'posn', PlotPosS(2)
case 'posu', PlotPosS(3)
case 'posm', PlotPosRcvM
case 'posr', PlotPosRcvR
case 'posd', PlotPosDate
case 'conv', PlotPosConv
case 'mapeh',PlotPosErrMap(1)
case 'mapev',PlotPosErrMap(0)
case 'maph', PlotPosVarMap(1)
case 'mapv', PlotPosVarMap(0)
case 'trpp', PlotTrp
end

% plot station position error --------------------------------------------------
function PlotPosH
data=get(gcf','userdata'); prm=data.prm;
clf, if isempty(data.err), return, end
i=find(strcmp(prm.rcv,prm.rcvs));
[t,j,k]=intersect(prm.time,data.timef);
err=data.err{i};
sig=data.sig{i};
posr=meann(data.posr{i});
if prm.est, err=err(j,:); sig=sig(j,:); end
axis off
ggt('newaxes','errv',[0.84,0.09,0.07,0.85],'move',[0,1,0,1],'xlim',[-1,1],...
    'ylim',prm.range,'yaxislocation','right','xtick',0,'xticklabel',{})
ylabel('Up (m)')
plot([-1,1],[0,0],'k'), plot(0,0,'k.','markersize',10);
if prm.ptype==1|prm.ptype==3, plot(zeros(size(err,1),1),err(:,7),'.','markersize',prm.psize(2)); end
if prm.ptype==2|prm.ptype==3, plot(zeros(size(err,1),1),err(:,7),'-','linewidth',prm.psize(1)); end
if prm.ptype==4, plot(zeros(size(err,1),1),err(:,7),'o','markersize',prm.psize(2)); end
if prm.showf(3)
    for n=1:length(data.timef)
        ep=mjdtocal(prm.td,data.timef(n));
        ggt('stext','errv',[0,err(n,7)],sprintf('%d/%d',ep(2:3)),1);
    end
end
ggt('newaxes','errh',[0.08,0.09,0.75,0.85],'equal','move',[1,1,1,1],...
    'xlim',prm.range*1.25,'ylim',prm.range)
xlabel('East (m)'), ylabel('North (m)')
plot([-100,100],[0,0],'k'), plot([0,0],[-100,100],'k')
plot(0,0,'k.','markersize',10);
if prm.ptype==1|prm.ptype==3, plot(err(:,5),err(:,6),'.','markersize',prm.psize(2)); end
if prm.ptype==2|prm.ptype==3, plot(err(:,5),err(:,6),'-','linewidth',prm.psize(1)); end
if prm.ptype==4, plot(err(:,5),err(:,6),'o','markersize',prm.psize(2)); end
if prm.showf(1)
    s=sprintf('REF: %s\nX: %13.4fm\nY: %13.4fm\nZ: %13.4fm',RefPos(prm.ref),posr);
    ggt('mtext','errh',s,1);
end
if prm.showf(2)
    f='MEAN E:%7.4fm N:%7.4fm U:%7.4fm\nRMS  E:%7.4fm N:%7.4fm U:%7.4fm';
    s=sprintf(f,meann(err(:,5:7)),rmsn(err(:,5:7)));
    ggt('mtext','errh',s,3);
end
if prm.showf(3)
    for n=1:length(data.timef)
        ep=mjdtocal(prm.td,data.timef(n));
        ggt('stext','errh',err(n,5:6),sprintf('%d/%d',ep(2:3)),3);
    end
end
ti=[tstr(prm.td,prm.time(1)),'-',tstr(prm.td,prm.time(end))];
title(['Receiver Position ',RcvName(prm.rcv,prm.showf(4)),' : ',ti]);

% plot station position 3d -----------------------------------------------------
function PlotPos3D
data=get(gcf','userdata'); prm=data.prm;
clf, if isempty(data.err), return, end
[fn,fs]=gut('getfont','g');
i=find(strcmp(prm.rcv,prm.rcvs));
err=data.err{i};
sig=data.sig{i};
axis off
tick=[prm.range(1),0,prm.range(2)];
ggt('newaxes','errh',[0.05,0.05,0.9,0.9],'equal','xlim',prm.range,'ylim',prm.range,...
    'zlim',prm.range,'xtick',tick,'ytick',tick,'ztick',tick)
ylabel('(m)')
if prm.ptype==1|prm.ptype==3
    plot3d(err(:,5),err(:,6),err(:,7),prm,'.','markersize',prm.psize(2));
end
if prm.ptype==2|prm.ptype==3
    plot3d(err(:,5),err(:,6),err(:,7),prm,'-','linewidth',prm.psize(1));
end
if prm.ptype==4
    plot3d(err(:,5),err(:,6),err(:,7),prm,'o','markersize',prm.psize(2));
end
if prm.showf(3)
    for n=1:length(prm.time)
        if mod(prm.time(n),60)==0
            plot3d(err(n,5),err(n,6),err(n,7),prm,'.','markersize',prm.psize(2)*2);
            ep=mjdtocal(prm.td,prm.time(n));
            text(err(n,5),err(n,6),err(n,7),sprintf('%d:%02d',ep(4:5)),...
                 'vertical','top','horizontal','center','fontname',fn,'fontsize',fs);
        end
    end
end
plot3(prm.range,[0,0],[0,0],'k')
plot3([0,0],prm.range,[0,0],'k')
plot3([0,0],[0,0],prm.range,'k')
plot3(0,0,0,'k.','markersize',prm.psize(2)*2);
text(prm.range(2)*1.1,0,0,'E','horizontal','center','fontname',fn,'fontsize',fs);
text(prm.range(1)*1.1,0,0,'W','horizontal','center','fontname',fn,'fontsize',fs);
text(0,prm.range(2)*1.1,0,'N','horizontal','center','fontname',fn,'fontsize',fs);
text(0,prm.range(1)*1.1,0,'S','horizontal','center','fontname',fn,'fontsize',fs);
ti=[tstr(prm.td,prm.time(1)),'-',tstr(prm.td,prm.time(end))];
title(['Receiver Position ',RcvName(prm.rcv,prm.showf(4)),' : ',ti]);
view(25,60), camproj('perspective'), rotate3d on

function plot3d(x,y,z,prm,varargin)
plot3(x,y,repmat(prm.range(1),size(x,1),1),varargin{:},'color',[0.8,0.8,0.8]);
%plot3(x,repmat(prm.range(2),size(x,1),1),z,varargin{:},'color',[0.8,0.8,0.8]);
%plot3(repmat(prm.range(1),size(x,1),1),y,z,varargin{:},'color',[0.8,0.8,0.8]);
plot3(x,y,z,varargin{:});

% plot station position error xyz ----------------------------------------------
function PlotPosX
data=get(gcf','userdata'); prm=data.prm;
clf, if isempty(prm.time)|isempty(data.err), return, end
i=find(strcmp(prm.rcv,prm.rcvs));
[t,j,k]=intersect(prm.time,data.timef);
time=prm.time;
err=data.err{i};
sig=data.sig{i};
posr=meann(data.posr{i});
if prm.est, time=time(j); err=err(j,:); sig=sig(j,:); end
label={'X (m)','Y (m)','Z (m)'};
margin=[0.08,0.03,0.04,0.012];
if length(prm.time)>1, xl=[prm.time(1),prm.time(end)]/3600;
else xl=[prm.time(1)-1,prm.time(1)+1]; end
for n=1:3
    if n<3, topts='nolabel'; else topts=''; end
    ggt('subplotv','err',n,3,'move',[1,1,1,1],'link',[1,0,1,0],...
        'taxis',prm.td,'topts',topts,'margin',margin,...
        'xlim',xl,'ylim',prm.range);
    plotval(time,err(:,n+1),prm.ptype,prm.psize);
    plotsig(time,sig(:,n+1),prm.nsig,prm.psize);
    s1=''; s2='';
    if prm.showf(1)
        s1=sprintf('REF: %s  %13.4fm ',RefPos(prm.ref),posr(n));
    end
    if prm.showf(2)
        s2=sprintf('MEAN: %6.4fm RMS: %6.4fm',meann(err(:,n+1)),rmsn(err(:,n+1)));
    end
    ggt('mtext',['err_',num2str(n)],[s1,s2],1); ylabel(label{n})
    if n==1
        ti=[tstr(prm.td,prm.time(1)),'-',tstr(prm.td,prm.time(end))];
        title(['Receiver Position ',RcvName(prm.rcv,prm.showf(4)),' : ',ti]);
    end
end

% plot station position enu ----------------------------------------------------
function PlotPosE
data=get(gcf','userdata'); prm=data.prm;
clf, if isempty(prm.time)|isempty(data.err), return, end
i=find(strcmp(prm.rcv,prm.rcvs));
[t,j,k]=intersect(prm.time,data.timef);
time=prm.time;
err=data.err{i};
sig=data.sig{i};
posr=meann(data.posr{i});
if prm.est, time=time(j); err=err(j,:); sig=sig(j,:); end
label={'East (m)','North (m)','Up (m)'};
margin=[0.08,0.03,0.04,0.012];
if length(prm.time)>1, xl=[prm.time(1),prm.time(end)]/3600;
else xl=[prm.time(1)-1,prm.time(1)+1]; end
for n=1:3
    if n<3, topts='nolabel'; else topts=''; end
    ggt('subplotv','err',n,3,'move',[1,1,1,1],'link',[1,0,1,0],...
        'taxis',prm.td,'topts',topts,'margin',margin,...
        'xlim',xl,'ylim',prm.range);
    plotval(time,err(:,n+4),prm.ptype,prm.psize);
    plotsig(time,sig(:,n+4),prm.nsig,prm.psize);
    s1=''; s2='';
    if prm.showf(1)
        s1=sprintf('REF: %s  ',RefPos(prm.ref));
    end
    if prm.showf(2)
        s2=sprintf('MEAN: %6.4fm RMS: %6.4fm',meann(err(:,n+4)),rmsn(err(:,n+4)));
    end
    ggt('mtext',['err_',num2str(n)],[s1,s2],1); ylabel(label{n})
    if n==1
        ti=[tstr(prm.td,prm.time(1)),'-',tstr(prm.td,prm.time(end))];
        title(['Receiver Position ',RcvName(prm.rcv,prm.showf(4)),' : ',ti]);
    end
end

% plot station position e/n/u --------------------------------------------------
function PlotPosS(n)
data=get(gcf','userdata'); prm=data.prm;
clf, if isempty(prm.time)|isempty(data.err), return, end
time=prm.time;
epos=readpos(prm.td,0,prm.rcvs,'','approx');
for m=1:length(prm.rcvs), gpos(m,:)=eceftogeod(epos(:,:,m)'); end
[lat,i]=sort(gpos(:,1));
[t,j,k]=intersect(prm.time,data.timef);
if prm.est, time=time(j); err=err(j,:); sig=sig(j,:); end
label={'East (m)','North (m)','Up (m)'};
margin=[0.11,0.03,0.04,0.012];
range=prm.range(2)-prm.range(1);
if length(prm.time)>1, xl=[prm.time(1),prm.time(end)]/3600;
else xl=[prm.time(1)-1,prm.time(1)+1]; end
ggt('subplotv','err',1,1,'move',[1,0,1,0],'taxis',prm.td,'margin',margin,...
    'xlim',xl,'ylim',[0,range*(length(prm.rcvs)+1)]);
err=[]; ytick=[]; ytlabel={};
for m=1:length(prm.rcvs)
    y=range*m; err=[err,y+data.err{i(m)}(:,n+4)];
    ytick=[ytick;y];
    ytlabel={ytlabel{:},RcvName(prm.rcvs{i(m)},prm.showf(4))};
end
plotval(time,err,prm.ptype,prm.psize);
set(gca,'xgrid','off','ytick',ytick,'yticklabel',ytlabel);
s1=''; if prm.showf(1), s1=sprintf('REF: %s  ',RefPos(prm.ref)); end
ggt('mtext','err_1',[s1,' GRID: ',num2str(range),'m'],1);

ti=[tstr(prm.td,prm.time(1)),'-',tstr(prm.td,prm.time(end))];
title(['Receiver Position Variations ',label{n},' : ',ti]);

% plot station position wave ---------------------------------------------------
function PlotPosW
data=get(gcf','userdata'); prm=data.prm;
clf, if isempty(prm.time)|isempty(data.err), return, end
[fn,fs]=gut('getfont','g');
i=find(strcmp(prm.rcv,prm.rcvs));
[t,j,k]=intersect(prm.time,data.timef);
time=prm.time;
err=data.err{i};
posr=meann(data.posr{i});
if prm.est, time=time(j); err=err(j,:); sig=sig(j,:); end
margin=[0.08,0.03,0.04,0.012];
range=prm.range(2)-prm.range(1);
if length(prm.time)>1, xl=[prm.time(1),prm.time(end)]/3600;
else xl=[prm.time(1)-1,prm.time(1)+1]; end
ggt('subplotv','err',1,1,'move',[1,0,1,0],'taxis',prm.td,'margin',margin,...
    'xlim',xl,'ylim',[-range*165,range*165]);
grid off;
label={'Up','North','East'};
ytick=[]; ytlabel={}; s1=''; s2='';
for n=1:3
    y=range*(n-2)*100;
    ytick=[ytick,y-range/2*100,y];
    plotval(time,err(:,8-n)*100+y,prm.ptype,prm.psize);
    text(prm.time(1)/3600,y+range*35,['  ',label{n}],'fontname',fn,'fontsize',fs,...
         'horizontal','left')
end
ytick=[ytick,y+range/2*100];
set(gca,'ytick',ytick);
ylabel('Displacements (cm)')
if prm.showf(1)
    s1=sprintf('REF: %s  X:%12.3fm Y:%12.3fm Z:%12.3fm \n',RefPos(prm.ref),posr);
end
if prm.showf(2)
    s2=sprintf('RMS: E %6.4fm, N %6.4fm U %6.4fm',rmsn(err(:,5:7)));
end
ggt('mtext','err_1',[s1,s2],1);
ti=[tstr(prm.td,prm.time(1)),'-',tstr(prm.td,prm.time(end))];
title(['Receiver Position ',RcvName(prm.rcv,prm.showf(4)),' : ',ti]);

% plot estimation mean error by stations ----------------------------------------
function PlotPosRcvM
data=get(gcf,'userdata'); prm=data.prm;
clf, if isempty(data.err), return, end
pos=[0.08,0.12,0.89,0.81];
[t,j,k]=intersect(prm.time,data.timef);
for n=1:length(prm.rcvs)
    e=data.err{n};
    if prm.est, err(n,:)=meann(e(j,5:7)); else err(n,:)=meann(e(:,5:7)); end
end
err=[err;meann(err)];
h=ggt('barplot',err,{prm.rcvs{:},'Average'},'ylim',[prm.range(1),prm.range(2)],'position',pos);
legend(h,{'East','North','Up'})
ylabel('Position Mean Error (m)');
ti=[tstr(prm.td,prm.time(1)),'-',tstr(prm.td,prm.time(end))];
title(['Receiver Position : ',ti,' (REF: ',RefPos(prm.ref),')']);

% plot estimation rms error by stations -----------------------------------------
function PlotPosRcvR
data=get(gcf,'userdata'); prm=data.prm;
clf, if isempty(data.err), return, end
pos=[0.08,0.12,0.89,0.81];
[t,j,k]=intersect(prm.time,data.timef);
for n=1:length(prm.rcvs)
    e=data.err{n};
    if prm.est, err(n,:)=rmsn(e(j,5:7)); else err(n,:)=rmsn(e(:,5:7)); end
end
err=[err;meann(err)];
h=ggt('barplot',err,{prm.rcvs{:},'Average'},'ylim',[0,prm.range(2)],'position',pos);
legend(h,{'East','North','Up'})
ylabel('Position RMS Error (m)');
ti=[tstr(prm.td,prm.time(1)),'-',tstr(prm.td,prm.time(end))];
title(['Receiver Position : ',ti,' (REF: ',RefPos(prm.ref),')']);

% plot estimatin error by date -------------------------------------------------
function PlotPosDate
data=get(gcf,'userdata'); prm=data.prm;
clf, if isempty(data.err), return, end
pos=[0.08,0.12,0.89,0.81];
[t,i,j]=intersect(prm.time,data.timef);
if prm.est, time=prm.time(i); else time=prm.time; end
err=[]; labels={};
for t=(floor(prm.time(1)/86400):floor(prm.time(end)/86400))*86400
    i=find(t<=time&time<t+86400); e=[];
    for n=1:length(prm.rcvs), e=[e;data.err{n}(i,5:7)]; end
    err=[err;rmsn(e)];
    d=mjdtocal(prm.td,t); labels={labels{:},sprintf('%d/%d',d(2:3))};
end
err=[err;meann(err)];
h=ggt('barplot',err,{labels{:},'Average'},'ylim',[0,prm.range(2)],'position',pos);
legend(h,{'East','North','Up'})
ylabel('Position RMS Error (m)');
ti=[tstr(prm.td,prm.time(1)),'-',tstr(prm.td,prm.time(end))];
title(['Receiver Position : ',ti,' (REF: ',RefPos(prm.ref),')']);

% summary station position estimation error ------------------------------------
function SummaryPosErr
data=get(gcf','userdata'); prm=data.prm;
if isempty(data.err), return, end
s={};
s=adds(s,'                RECEIVER POSITIONS : POSITION ERROR (wrt %s)',RefPos(prm.ref));
s=adds(s,'');
s=adds(s,' TIME    : %s-%s GPST',tstr(prm.td,prm.time(1)),tstr(prm.td,prm.time(end)));
s=adds(s,'                mean error (m)        rms  error (m)       repeatability (m)');
s=adds(s,'            east   north     up     east   north   up     east   north   up');
s=adds(s,' -----------------------------------------------------------------------------');
if prm.est, [t,i,j]=intersect(prm.time,data.timef); else i=1:length(prm.time); end
for n=1:length(prm.rcvs)
    if ~isempty(data.err{n})
        e(n,:)=[meann(data.err{n}(i,5:7)),rmsn(data.err{n}(i,5:7)),stdn(data.err{n}(i,5:7))];
    else
        e(n,:)=repmat(nan,1,9);
    end
    rcv=prm.rcvs{n};
    if prm.showf(4)==2, s1=[' : ',RcvName(rcv,prm.showf(4))]; else s1=''; end
    if prm.adj&any(strcmp(rcv,data.rcvrf)), rcv=[rcv,'(R)']; end
    s=adds(s,' %-7s :%7.4f %7.4f %7.4f  %6.4f %6.4f %6.4f  %6.4f %6.4f %6.4f %s',rcv,e(n,:),s1);
end
s=adds(s,' -----------------------------------------------------------------------------');
s=adds(s,' average :%7.4f %7.4f %7.4f  %6.4f %6.4f %6.4f  %6.4f %6.4f %6.4f',meann(e));
s=adds(s,' median  :%7.4f %7.4f %7.4f  %6.4f %6.4f %6.4f  %6.4f %6.4f %6.4f',mediann(e));
if ~isempty(prm.adj)
    s=adds(s,'                                                         (R):%s Station',RefPos(prm.adj));
end
if ~isempty(data.tr)
    s=adds(s,'');
    s=adds(s,'                     COORDINATE TRANSFORMATION PARAMETERS');
    s=adds(s,'              dx(m)    dy(m)    dz(m)   Rx(mas)  Ry(mas)  Rz(mas) scale(ppb)');
    s=adds(s,' -----------------------------------------------------------------------------');
    for n=1:size(data.tr,1)
        t=caltomjd(data.tr(n,1:3));
        if prm.td+prm.time(1)/86400<=t&t<=prm.td+prm.time(end)/86400
            s=adds(s,' %04d/%02d/%02d:%8.4f %8.4f %8.4f %8.4f %8.4f %8.4f %8.4f',data.tr(n,:));
        end
    end
end
ti=prm.dirs.est;
Viewer(s,['Summary of Receiver Position : ',ti]);

% plot transformation parameters -----------------------------------------------
function PlotTrp
data=get(gcf','userdata'); prm=data.prm;
clf, if isempty(data.tr), return, end
for n=1:size(data.tr,1), t(n)=(caltomjd(data.tr(n,1:3))-prm.td)*24; end
label={'dx/dy/dz (m)','Rx/Ry/Rz (mas)','Scale (ppb)'};
ls={{'dx','dy','dz'},{'Rx','Ry','Rz'}};
margin=[0.08,0.03,0.04,0.012];
for n=1:3
    if n<3, topts='nolabel'; tr=data.tr(:,n*3+1:n*3+3);
    else topts=''; tr=data.tr(:,end); end
    ggt('subplotv','err',n,3,'move',[1,1,1,1],'link',[1,0,1,0],...
        'taxis',prm.td,'topts',topts,'margin',margin,...
        'xlim',[prm.time(1),prm.time(end)]/3600);
    plot(t,tr,'-','linewidth',prm.psize(1));
    plot(t,tr,'.','markersize',prm.psize(2));
    ylabel(label{n});
end
ti=[tstr(prm.td,prm.time(1)),'-',tstr(prm.td,prm.time(end))];
gut('setstring','title',['Transformation Parameters : ',ti]);

% output transformation parameters --------------------------------------------
function OutputTrp
data=get(gcf','userdata'); prm=data.prm;
if isempty(data.pos)|isempty(prm.adj), return, end
s={};
s=adds(s,'%%              TRANSFORMATION PARAMETERS ADJUSTING TO %s',RefPos(prm.adj));
s=adds(s,'%%             dx(m)    dy(m)    dz(m)   Rx(mas)  Ry(mas)  Rz(mas) scale(ppb)');
s=adds(s,'%%-----------------------------------------------------------------------------');
s=adds(s,'function prm=prm_transprm');
s=adds(s,'prm=[');
for n=1:size(data.tr,1)
    s=adds(s,'%4d,%2d,%2d,%8.4f,%8.4f,%8.4f,%8.4f,%8.4f,%8.4f,%8.4f',data.tr(n,:));
end
s=adds(s,'];');
ti=prm.dirs.est;
Viewer(s,['Transformation Parameters : ',ti]);

% plot station position error map ---------------------------------------------
function PlotPosErrMap(hv)
data=get(gcf','userdata'); prm=data.prm;
clf, if isempty(data.err), return, end
[fname,fsize]=gut('getfont','g');
p=get(gcf,'position'); pos=[0.001,0.001,0.998,(p(4)-20)/p(4)];
axis off, box on
gmt('mmap','proj',prm.map.proj,'cent',prm.map.cent,'base',prm.map.base,...
    'scale',prm.map.scale,'pos',pos,'fontname',fname,'fontsize',fsize);
gmt('mcoast','lcolor',prm.map.color{1},'scolor',prm.map.color{2},'ccolor',prm.map.color{3});
gmt('mgrid','gint',prm.map.gint(1),'lint',prm.map.gint(2),'color',prm.map.color{4});
[xl,yl]=gmt('getlim');
if prm.showf(1)
    msg=['REF: ',RefPos(prm.ref)];
    text(xl(2)*0.95,yl(2)*0.98,msg,'vertical','top','horizontal','right',...
         'fontname',fname,'fontsize',fsize);
end
range=(prm.range(2)-prm.range(1))/10;
label=[num2str(range),'m'];
if hv
    type='Horizontal'; color='b';
    [lon1,lat1]=gmt('xytoll',xl(2)*0.8,yl(1)*0.85);
    [lon2,lat2]=gmt('xytoll',xl(2)*0.85,yl(1)*0.85);
    gmt('marrow',lon1,lat1,5,0,'color',color,'linewidth',prm.psize(1));
    gmt('mtext',lon2,lat2,label,'vertical','bottom','horizontal','center','fontname',fname,'fontsize',fsize);
else
    type='Vertical'; color='r';
    [lon1,lat1]=gmt('xytoll',xl(2)*0.8,yl(1)*0.92);
    [lon2,lat2]=gmt('xytoll',xl(2)*0.8,yl(1)*0.85);
    gmt('marrow',lon1,lat1,0,5,'color',color,'linewidth',prm.psize(1));
    h=gmt('mtext',lon2,lat2,[' ',label],'vertical','middle','horizontal','left','fontname',fname,'fontsize',fsize);
end
if prm.est, [t,i,j]=intersect(prm.time,data.timef); else i=1:length(prm.time); end
for n=1:length(prm.rcvs)
    if ~isempty(data.pos{n})
        gpos=eceftogeod(meann(data.pos{n}(i,:))');
        gmt('mplot',gpos(2),gpos(1),color,'marker','.','markersize',prm.psize(2));
        if prm.showf(4)>0
            rcv=RcvName(prm.rcvs{n},prm.showf(4));
            h=gmt('mtext',gpos(2),gpos(1),rcv,'horizontal','center','vertical','top');
            if isstruct(prm.rfont), set(h,prm.rfont);
            else set(h,'fontname',fname,'fontsize',fsize); end
        end
        err=meann(data.err{n}(i,:));
        if hv, gmt('marrow',gpos(2),gpos(1),err(5)*5/range,err(6)*5/range,'b','linewidth',prm.psize(1));
        else gmt('marrow',gpos(2),gpos(1),0,err(7)*5/range,'color','r','linewidth',prm.psize(1)); end
    end
end
ti=[tstr(prm.td,prm.time(1)),'-',tstr(prm.td,prm.time(end))];
title(['Receiver Position ',type,' Error : ',ti]);

% plot station position variation map ------------------------------------------
function PlotPosVarMap(hv)
data=get(gcf','userdata'); prm=data.prm;
clf, if isempty(data.err), return, end
[fname,fsize]=gut('getfont','g');
p=get(gcf,'position'); pos=[0.001,0.001,0.998,(p(4)-20)/p(4)];
axis off, box on
gmt('mmap','proj',prm.map.proj,'cent',prm.map.cent,'base',prm.map.base,...
    'scale',prm.map.scale,'pos',pos,'fontname',fname,'fontsize',fsize);
gmt('mcoast','lcolor',prm.map.color{1},'scolor',prm.map.color{2},'ccolor',prm.map.color{3});
gmt('mgrid','gint',prm.map.gint(1),'lint',prm.map.gint(2),'color',prm.map.color{4});
[xl,yl]=gmt('getlim');
range=(prm.range(2)-prm.range(1))/10;
label=[num2str(range),'m'];
if hv
    type='Horizontal'; color='b';
    [lon1,lat1]=gmt('xytoll',xl(2)*0.8,yl(1)*0.85);
    [lon2,lat2]=gmt('xytoll',xl(2)*0.85,yl(1)*0.85);
    gmt('marrow',lon1,lat1,5,0,'color',color,'linewidth',prm.psize(1));
    gmt('mtext',lon2,lat2,label,'vertical','bottom','horizontal','center','fontname',fname,'fontsize',fsize);
else
    type='Vertical'; color='r';
    [lon1,lat1]=gmt('xytoll',xl(2)*0.8,yl(1)*0.92);
    [lon2,lat2]=gmt('xytoll',xl(2)*0.8,yl(1)*0.85);
    gmt('marrow',lon1,lat1,0,5,'color',color,'linewidth',prm.psize(1));
    gmt('mtext',lon1,lat2,label,'vertical','middle','horizontal','left','fontname',fname,'fontsize',fsize);
end
if prm.est, [t,i,j]=intersect(prm.time,data.timef); else i=1:length(prm.time); end
for n=1:length(prm.rcvs)
    p=data.pos{n}(end,:);
    gpos=eceftogeod(data.pos{n}(i(1),:)');
    gmt('mplot',gpos(2),gpos(1),color,'marker','.');
    if prm.showf(4)>0
        rcv=RcvName(prm.rcvs{n},prm.showf(4));
        h=gmt('mtext',gpos(2),gpos(1),rcv,'horizontal','center','vertical','top');
        if isstruct(prm.rfont), set(h,prm.rfont);
        else set(h,'fontname',fname,'fontsize',fsize); end
    end
    [p,E]=geodtoecef(gpos);
    dp=(data.pos{n}(i(end),:)-data.pos{n}(i(1),:))*E';
    if hv, gmt('marrow',gpos(2),gpos(1),dp(1)*5/range,dp(2)*5/range,'linewidth',prm.psize(1));
    else gmt('marrow',gpos(2),gpos(1),0,dp(3)*5/range,'color','r','linewidth',prm.psize(1)); end
end
ti=[tstr(prm.td,prm.time(1)),' -> ',tstr(prm.td,prm.time(end))];
title(['Receiver Position ',type,' Variation : ',ti]);

% read data --------------------------------------------------------------------
function ReadData
data=get(gcf','userdata'); prm=data.prm;
[prm.td,ts]=caltomjd(prm.tstart); [tn,te]=caltomjd(prm.tend);
time=ts:prm.tint:te+(tn-prm.td)*86400; if isempty(time), return, end
prm.time=time;
data.pos={}; data.err={}; data.sig={}; var={};
gut('newmsgbar','msgbar','',[480,50],[0,0.7,1]);
gut('setmsgbar','msgbar','reading receiver positions',0);
[poss,covs]=readpos(prm.td,prm.time,prm.rcvs,prm.dirs.est,prm.fb,prm.tunit);

if prm.srday(1)<=prm.srday(2)&prm.srday(2)<0 % sidereal filter
    for n=1:length(prm.rcvs)
        poss(:,:,n)=srfilt(prm.td,prm.time,prm.rcvs{n},poss(:,:,n),prm.srday,...
                           prm.srprm,prm.dirs.est,prm.fb,prm.tunit);
    end
end
for n=1:length(prm.rcvs)
    data.pos{n}=poss(:,:,n); var{n}=covs(:,:,n);
end
ti=['Receiver Position : ',prm.dirs.est,' (',prm.fb,')'];
tu=prm.tunit*3600;
ts=(floor(prm.time(1)/tu):floor(prm.time(end)/tu))*tu;

data.timef=[];
for n=1:length(ts)
    i=find(ts(n)<=prm.time&prm.time<ts(n)+prm.tunit*3600);
    switch prm.fb
    case 'posf', j=i(end);
    case 'posb', j=i(1);
    otherwise, j=floor((i(end)+i(1))/2);
    end
    data.timef(n,1)=prm.time(j);
    for m=1:length(prm.rcvs), posf(m,:,n)=data.pos{m}(j,:); end
end
% helmert transformation parameters adjusting to reference frame
data.tr=[]; data.rcvrf={};
if ~isempty(prm.rcvs)
    if ~isempty(prm.adj)
        for n=1:length(ts)
            % rf positions
            posrf=shiftdim(readpos(prm.td,data.timef(n),prm.rcvs,prm.dirs.pos,prm.adj),1)';
            data.rcvrf=prm.rcvs(~isnan(posrf(:,1)));
            
            % transformation parameters (est. positions->rf positions)
            ep=mjdtocal(prm.td,ts(n));
            data.tr(n,:)=[ep(1:3),esthelmert(posf(:,:,n),posrf,4)'];
        end
    elseif ~isempty(prm.trprm)
        [d,f,e]=fileparts(prm.trprm);
        if strcmp(e,'.m'), wd=pwd; cd(d); data.tr=feval(f); cd(wd); end
    end
end
% transformation corrections
if ~isempty(data.tr)
    for n=1:size(data.tr,1)
        t=(caltomjd(data.tr(n,1:3))-prm.td)*86400;
        i=find(t<=prm.time&prm.time<t+86400);
        for m=1:length(prm.rcvs)
            data.pos{m}(i,:)=helmert(data.pos{m}(i,:),data.tr(n,4:end)');
        end
    end
end
gut('setmsgbar','msgbar','reading reference positions',[]);

% read references and compute errors
data.posr=ReadRef(prm.td,prm.time,prm.rcvs,prm.dirs,prm.ref,data.pos,prm.fb,prm.tunit);
[data.err,data.sig]=ErrEstPos(prm.td,prm.time,prm.rcvs,prm.dirs,data.pos,var,...
                              data.posr,prm.tunit,prm.bpfilt);
gut('closemsgbar','msgbar');
data.prm=prm;
set(gcf,'name',ti,'userdata',data);

% update station menu
delete(gut('geth','rcvs'));
gut('newmenu','rcvs','&Receivers',prm.rcvs,prm.rcvs,[mfilename,' cbRcv']);
gut('chkmenu','rcvs',prm.rcv);

% read reference positions -----------------------------------------------------
function posr=ReadRef(td,time,rcvs,dirs,ref,pose,fb,tunit)
posr={};
switch ref
case {'estm','ests','este','est0','estf','estp'}
    for n=1:length(rcvs)
        switch ref
        case 'estm', posr{n}=repmat(meann(pose{n}),length(time),1);
        case 'ests', posr{n}=repmat(pose{n}(1,:),length(time),1);
        case 'este', posr{n}=repmat(pose{n}(end,:),length(time),1);
        case 'est0'
            tu=tunit*3600; tint=time(2)-time(1);
            t=floor(time(1)/tu)*tu:tint:ceil(time(end)/tu)*tu-tint;
            pos0=readpos(td,t,rcvs{n},dirs.est,fb,tunit);
            posr{n}=repmat(mean(pos0(~isnan(pos0(:,1)),:)),length(time),1);
        case 'estf' % linear fit
            posr{n}=repmat(nan,length(time),3);
            for m=1:3
                if size(pose{n},1)>1
                    i=find(~isnan(pose{n}(:,m)));
                    x=[time(i)',ones(length(i),1)]\pose{n}(i,m);
                    posr{n}(i,m)=x(1)*time(i)'+x(2);
                else
                    posr{n}(:,m)=pose{n}(:,m);
                end
            end
        case 'estp' % poly fit
            posr{n}=repmat(nan,length(time),3);
            for m=1:3
                i=find(~isnan(pose{n}(:,m)));
                if length(i)>2
                    [p,s,mu]=polyfit(time(i)',pose{n}(i,m),2);
                    posr{n}(i,m)=polyval(p,time(i)',[],mu);
                end
            end
        end
    end
otherwise
    pos=readpos(td,time,rcvs,dirs.pos,ref);
    for n=1:length(rcvs), posr{n}=pos(:,:,n); end
end

% sidereal filter --------------------------------------------------------------
function poss=siderealfilt(poss,prm)
if prm.srday(2)<prm.srday(1), return, end
srsec=86164-prm.srday(3); % sidereal day (sec)
day=prm.srday(1):prm.srday(2);
off=zeros(length(prm.time),3,length(prm.rcvs)); nn=zeros(length(prm.rcvs));
for n=1:length(day)
    time=prm.time+round(srsec*day(n)/prm.tint)*prm.tint;
    pos=readpos(prm.td,time,prm.rcvs,prm.dirs.est,prm.fb,prm.tunit);
    for m=1:length(prm.rcvs)
        corr=pos(:,:,m);
        i=find(~isnan(corr(:,1)));
        if ~isempty(i)
            % delete dc
            corr(i,:)=detrend(corr(i,:),'constant');
            
            % interpolate/extrapolate outage period
            if isnan(corr(1,1)), corr(1,:)=corr(i(1),:); i=[1;i]; end
            if isnan(corr(end,1)), corr(end,:)=corr(i(end),:); i=[i;length(time)]; end
            corr=interp1(time(i),corr(i,:),time);
            
            % bandbass filter
            if prm.srprm(1)<prm.srprm(2)
                corr=bpf(time,corr,prm.srprm);
            end
            nn(m)=nn(m)+1; off(:,:,m)=(off(:,:,m)*(nn(m)-1)+corr)/nn(m);
        end
    end
end
poss=poss-off;

% position estimation errors ---------------------------------------------------
function [err,sig]=ErrEstPos(td,time,rcvs,dirs,pos,var,posr,tunit,filt)
err={}; sig={};
for n=1:length(rcvs)
    err{n}=zeros(length(time),7); sig{n}=err{n};
    for m=1:length(time)
        [p,E]=geodtoecef(eceftogeod(posr{n}(m,:)'));
        errxyz=pos{n}(m,:)-posr{n}(m,:);
        errenu=errxyz*E';
        varenu=diag(E*diag(var{n}(m,:))*E')';
        err{n}(m,:)=[norm(errenu(1:2)),errxyz,errenu];
        sig{n}(m,:)=sqrt([sum(varenu(1:2)),var{n}(m,:),varenu]);
    end
    if filt(1)<filt(2)
        err{n}=bpfilt(time,err{n},filt(1),filt(2));
    end
end

% station name -----------------------------------------------------------------
function str=RcvName(rcv,showf)
persistent prmr, if isempty(prmr), prmr=prm_gpsrcvs; end
if any(showf==[0,1,3]), str=[rcv,' ']; else str=''; end
if any(showf==[2,3])
    i=min(find(strcmp(rcv,prmr(:,1))));
    if ~isempty(i), str=[str,prmr{i,3},' ']; end
end

% reference point --------------------------------------------------------------
function str=RefPos(ref)
strs={'Average','Linear Fit','Poly Fit','Start Pos','End Pos','Est Average','ITRF97','ITRF2000','IGS Final','GSIPOS','IGS00','IGb00'};
refs={'estm','estf','estp','ests','este','est0','itrf97','itrf2000','igssnx','gsipos','igs00','igb00'};
i=find(strcmp(ref,refs)); if isempty(i), str=''; else str=strs{i}; end

% show viewer ------------------------------------------------------------------
function Viewer(str,ti)
h=gut('newviewer','',ti,[600,419,0,-48],'str',str);

% plot value -------------------------------------------------------------------
function h=plotval(t,x,ptype,psize)
if ptype==1|ptype==3, h=plot(t/3600,x,'.','markersize',psize(2)); end
if ptype==2|ptype==3, h=plot(t/3600,x,'-','linewidth',psize(1)); end
if ptype==4, h=plot(t/3600,x,'o','markersize',psize(2)); end

% plot sigma -------------------------------------------------------------------
function plotsig(t,sig,nsig,psize)
if nsig>0
    plot(t/3600,-nsig*sig,'r:','linewidth',psize(1))
    plot(t/3600, nsig*sig,'r:','linewidth',psize(1))
end

% mean/std without nan ---------------------------------------------------------
function m=meann(x)
for n=1:size(x,2)
    i=find(~isnan(x(:,n)));
    if isempty(i), m(1,n)=nan; else m(1,n)=mean(x(i,n),1); end
end
function s=mediann(x)
for n=1:size(x,2)
    i=find(~isnan(x(:,n)));
    if isempty(i), s(1,n)=nan; else s(1,n)=median(x(i,n),1); end
end
function s=rmsn(x)
for n=1:size(x,2)
    i=find(~isnan(x(:,n)));
    if isempty(i), s(1,n)=nan; else s(1,n)=sqrt(mean(x(i,n).^2,1)); end
end
function s=stdn(x)
for n=1:size(x,2)
    i=find(~isnan(x(:,n)));
    if isempty(i), s(1,n)=nan; else s(1,n)=std(x(i,n),1); end
end

% date/time string -------------------------------------------------------------
function str=tstr(td,ts)
str=sprintf('%04d/%02d/%02d %02d:%02d:%02.0f',mjdtocal(td,ts));

% add string -------------------------------------------------------------------
function s=adds(s,varargin), s={s{:},sprintf(varargin{:})};
