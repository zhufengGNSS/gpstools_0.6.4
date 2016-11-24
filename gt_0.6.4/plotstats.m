function plotstats(varargin)
%-------------------------------------------------------------------------------
% [system] : GpsTools
% [module] : plot residuals
% [func]   : plot residuals
% [argin]  : 'prm',prm  = parameters
% [argout] : h = figurea handle
% [note]   :
% [version]: $Revision: 12 $ $Date: 2008-11-25 10:02:15 +0900 (火, 25 11 2008) $
%            Copyright(c) 2004-2006 by T.Takasu, all rights reserved
% [history]: 04/11/17  0.1  new
%-------------------------------------------------------------------------------
if nargin>0&strncmp(varargin{1},'cb',2), feval(varargin{:});
else h=PlotMain(varargin{:}); end

% main -------------------------------------------------------------------------
function h=PlotMain(varargin)
h=gut('newfigg','','Residuals/Statistics',[600,400,0,-68],[mfilename,' cbClose']');
if isempty(h), return, end
prm=loadprm('prm_plotstats','prm_plotstats_def');
if nargin>=2&strcmp(varargin{1},'prm'), prm=varargin{2}; end
cb1=[mfilename,' cbPlot1'];
cb2=[mfilename,' cbPlot2'];
cb3=[mfilename,' cbPlot3'];
gut('newmenu','data','&Data ',...
    {'&Read...','&New Window','-','&Export Plot...','-','&Options...'},{},...
    {[mfilename,' cbRead'],[mfilename,' cbNew'],'',[mfilename,' cbExport'],'',...
     [mfilename,' cbOpt']});
gut('newmenu','plot','&Plot',...
    {'Residuals','Residuals East/North/Up','Residuals Az/El',...
     'Residuals Heading/Nadir','Residuals Sky Map','Residulas 3D Sky Map',...
     'Residuals Average Sky Map','Residuals Std Dev Sky Map','-','Statistics...',...
     'Sat Antenna PCV'},...
     {'res','enu','azel','thph','sky','sky3','skya','skys','','',''},...
     {cb1,cb1,cb1,cb1,cb1,cb1,cb1,cb1,'',cb2,cb3});
gut('newmenu','sats','&Satellite',{'ALL'},{'ALL'},[mfilename,' cbSat']);
gut('newmenu','rcvs','&Receiver',{'ALL'},{'ALL'},[mfilename,' cbRcv']);
gut('chkmenu','plot',prm.type);
gut('chkmenu','sats','ALL');
gut('chkmenu','rcvs','ALL');
set(gut('newbtnh','sbtn',[0,0,72,16],{'<','>','<','>'},[mfilename,' cbSwBtn']),...
    'fontsize',8);
set(gcf,'resizefcn',[mfilename,' cbResize']);
LocatePos;
data.prm=prm; data.t=[]; set(h,'userdata',data);
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
data=get(gcf,'userdata'); prm=data.prm;
sel1={{'Estimated(Forward)','Estimated(Backward)'},{'resf','resb'}};
sel2={{'Combined','RINEX','RINEX(3H)','RINEX(1H)','RINEX(15min)'},...
      {'brdc','rinex','rinex3','rinex1','rinexh'}};
prm1={
't','Start Time (GPST)',[]
't','End Time (GPST)',  []
'n','Processing Unit Time (hr)', 24
'b','Receivers',[mfilename,' cbSelRcv']
'p','Residuals/Statistics', sel1
'p','Navigation Messages', sel2
};
prm2={
' ','Estimation Data Directory',''
'd','',prm.dirs.est
' ','Navigation Message Directory',''
'd','',prm.dirs.nav
};
gut('newdlg','','Read Data',[314,268]);
gut('newprms','prm1',[12,114,295,25,118],prm1);
gut('newprms','prm2',[12,34,295,18],prm2);
gut('newokcancelbtn','',[127,4,180,23]);
gut('setprms','prm1',prm.tstart,prm.tend,prm.tunit,'',prm.fb,prm.srcnav);
gut('setudata','prm1_4',prm.rcvs);
if ~gut('waitok'), return, end
[prm.tstart,prm.tend,prm.tunit,q,prm.fb,prm.srcnav]=gut('getprms','prm1');
[q,prm.dirs.est,q,prm.dirs.nav]=gut('getprms','prm2');
prm.rcvs=gut('getudata','prm1_4');
close
data.prm=prm; set(gcf,'userdata',data);
SaveSetting; ReadData; UpdatePlot;

% callback on menu select receiver ---------------------------------------------
function cbSelRcv
[rcvs,ok]=editlist('Receivers',get(gcbo,'userdata'),'rcv');
if ok, set(gcbo,'userdata',rcvs); end

% callback on menu export plot -------------------------------------------------
function cbExport
persistent file dirs
if isempty(file), file='plotstats'; end, if isempty(dirs), dirs=''; end
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
if strcmp(prm.type,'res')
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
sel2={{'Line','Dot','Line and Dot','Circle'},1:4};
sel3={{'OFF','Stack','Correct'},0:2};
data=get(gcf,'userdata'); prm=data.prm;
prm1={
'p','Show Mean and RMS Residuals',sel1
'p','Plot Average of Residuals',  sel1
'p','Plot Std. Dev. of Residuals',sel1
'p','Plot All Residuals on Background',sel1
'p','Exclude Eclipsing Satellite',sel1
's','3D Skyplot View Angle (deg)',[45,30]
'p','Plot Style',                 sel2
'n','Plot Dot Size',              5
's','Mesh Interval Az/El (deg)',  [5,5]
's','Plot Axis Range',            [-0.1,0.1]
};
gut('newdlg','','Options',[320,288]);
gut('newprms','prm1',[15,34,300,25,120],prm1);
gut('setprms','prm1',prm.showf(1),prm.showf(2),prm.showf(3),prm.showf(4),...
    prm.execl,prm.vangle,prm.ptype,prm.psize,prm.mesh,prm.range);
gut('newokcancelbtn','',[135,4,180,23]);
if ~gut('waitok'), return, end
[prm.showf(1),prm.showf(2),prm.showf(3),prm.showf(4),prm.execl,prm.vangle,...
 prm.ptype,prm.psize,prm.mesh,prm.range]=gut('getprms','prm1');
close
data.prm=prm; set(gcf,'userdata',data)
SaveSetting; UpdatePlot;

% callback on menu plot --------------------------------------------------------
function cbPlot1
data=get(gcf,'userdata'); prm=data.prm;
gut('unchkmenu','plot');
prm.type=get(gcbo,'userdata');
gut('chkmenu','plot',prm.type);
data.prm=prm; set(gcf,'userdata',data)
UpdatePlot

function cbPlot2, SummaryStats
function cbPlot3, SatAntPcv

% callback on menu satellite ---------------------------------------------------
function cbSat
data=get(gcf,'userdata'); prm=data.prm;
prm.sat=get(gcbo,'userdata');
gut('unchkmenu','sats');
gut('chkmenu','sats',prm.sat);
data.prm=prm; set(gcf,'userdata',data);
UpdatePlot;

% callback on menu receiver ----------------------------------------------------
function cbRcv
data=get(gcf,'userdata'); prm=data.prm;
prm.rcv=get(gcbo,'userdata');
gut('unchkmenu','rcvs');
gut('chkmenu','rcvs',prm.rcv);
data.prm=prm; set(gcf,'userdata',data);
UpdatePlot;

% callback on push swithch button ----------------------------------------------
function cbSwBtn
data=get(gcf,'userdata'); prm=data.prm;
i=find(strcmp(prm.sats,prm.sat));
j=find(strcmp(prm.rcvs,prm.rcv));
switch get(gcbo,'tag')
case 'sbtn_1', if ~isempty(i), prm.sat=prm.sats{max(i-1,1)}; end
case 'sbtn_2', if ~isempty(i), prm.sat=prm.sats{min(i+1,length(prm.sats))}; end
case 'sbtn_3', if ~isempty(j), prm.rcv=prm.rcvs{max(j-1,1)}; end
case 'sbtn_4', if ~isempty(j), prm.rcv=prm.rcvs{min(j+1,length(prm.rcvs))}; end
end
gut('unchkmenu','sats'); gut('chkmenu','sats',prm.sat);
gut('unchkmenu','rcvs'); gut('chkmenu','rcvs',prm.rcv);
data.prm=prm; set(gcf,'userdata',data);
UpdatePlot;

% locate position --------------------------------------------------------------
function LocatePos
p=get(gcf,'position'); m=15;
gut('setpos','sbtn_1',[p(3)-m-64,p(4)-16]);
gut('setpos','sbtn_2',[p(3)-m-48,p(4)-16]);
gut('setpos','sbtn_3',[p(3)-m-32,p(4)-16]);
gut('setpos','sbtn_4',[p(3)-m-16,p(4)-16]);

% save settingsrs --------------------------------------------------------------
function SaveSetting
data=get(gcf,'userdata'); prm=data.prm;
[path,file]=fileparts(which(mfilename));
save(fullfile(path,'settings','prm_plotstats.mat'),'prm');

% update plot ------------------------------------------------------------------
function UpdatePlot
data=get(gcf,'userdata'); prm=data.prm;
switch prm.type
case 'res',  PlotResidual;
case 'enu',  PlotResidualEnu;
case 'azel', PlotResidualAzEl;
case 'thph', PlotResidualThPh;
case 'sky',  PlotResidualSky;
case 'sky3', PlotResidualSky3D;
case 'skya', PlotResidualSkyAve(0);
case 'skys', PlotResidualSkyAve(1);
otherwise,   PlotEstStat;
end

% plot residulas ---------------------------------------------------------------
function PlotResidual
data=get(gcf,'userdata'); prm=data.prm;
if isempty(data.t), return, end
time=data.t/3600;
index=data.index;
res=data.ress(:,2);
out=data.outs(:);
ecl=data.ecl;
if prm.execl
    i=find(ecl); time(i)=[]; index(i,:)=[]; res(i)=[]; out(i)=[]; m=[]; ecl=[];
end
i=ones(length(res),1);
if ~strcmp(prm.sat,'ALL')
    n=find(strcmp(prm.sat,prm.sats)); if isempty(n), return, end
    i=i&any(index(:,1:2)==n,2);
end
if ~strcmp(prm.rcv,'ALL')
    n=find(strcmp(prm.rcv,prm.rcvs)); if isempty(n), return, end
    i=i&any(index(:,3:4)==n,2);
end
j=find(i&~out); k=find(i&out); m=find(i&ecl);
clf
margin=[0.08,0.03,0.04,0.012];
ggt('subplotv','res',1,1,'move',[1,1,1,1],'link',[1,0,1,0],'taxis',prm.td,...
    'margin',margin,'xlim',[prm.time(1),prm.time(end)]/3600,'ylim',prm.range);
if prm.ptype==4, mark='o'; else mark='.'; end
if prm.showf(4)
    plot(time,res,'.','color',[0.8,0.8,0.8],'marker',mark,'markersize',prm.psize)
end
if prm.ptype>1
    plot(time(j),res(j,:),'b','linestyle','none','marker',mark,'markersize',prm.psize)
end
plot(time(k),res(k,:),'m','linestyle','none','marker',mark,'markersize',prm.psize)
if prm.ptype==1|prm.ptype==3
    plotsep(time(j),time(j),res(j),'b-',prm.tint);
    plotsep(time(m),time(m),res(m),'r-',prm.tint);
end
if prm.showf(1)
    s=sprintf('MEAN : %.4fm RMS : %.4fm',meann(res(j,:)),rmsn(res(j,:)));
    ggt('mtext','res_1',s,1);
end
if any(prm.showf(2:3))
    t=prm.time(1)/3600:prm.time(end)/3600+1; rsm=repmat(nan,length(t),1); rss=rsm;
    for n=1:length(t)
        i=find(t(n)-0.5<=time(j)&time(j)<t(n)+0.5);
        rsm(n)=meann(res(j(i)));
        rss(n)=stdn(res(j(i)));
    end
    if prm.showf(2), plot(t,rsm,'r-'), plot(t,rsm,'r.'), end
    if prm.showf(3), plot(t,rss,'m-'), plot(t,rss,'m.'), end
end
ylabel('Residuals (m)');
ti=sprintf('Postfit Residuals %s-%s : %s',prm.sat,prm.rcv,pstr(prm.td,prm.time));
title(ti);

% plot residulas e/n/u ----------------------------------------------------------
function PlotResidualEnu
data=get(gcf,'userdata'); prm=data.prm;
if isempty(data.t), return, end
if strcmp(prm.rcv,'ALL'), rcv=prm.rcvs{1}; else rcv=prm.rcv; end
n=find(strcmp(rcv,prm.rcvs)); if isempty(n), return, end
i=find(any(data.index(:,3:4)==n,2));
time=data.t(i)/3600;
index=data.index(i,:);
res=data.ress(i,2);
out=data.outs(i);
azels=data.azels(i,:);
i=ones(length(res),1);
if ~strcmp(prm.sat,'ALL')
    n=find(strcmp(prm.sat,prm.sats)); if isempty(n), return, end
    i=i&any(index(:,1:2)==n,2);
end
j=find(i&~out); k=find(i&out);
clf
res=[res.*cos(azels(:,2)).*sin(azels(:,1)),res.*cos(azels(:,2)).*cos(azels(:,1)),...
     res.*sin(azels(:,2))];
label={'East (m)','North (m)','Up (m)'};
margin=[0.08,0.03,0.04,0.012];
for n=1:3
    if n<3, topts='nolabel'; else topts=''; end
    ggt('subplotv','res',n,3,'move',[1,1,1,1],'link',[1,0,1,0],...
        'taxis',prm.td,'topts',topts,'margin',margin,...
        'xlim',[prm.time(1),prm.time(end)]/3600,'ylim',prm.range);
    if prm.ptype==4, mark='o'; else mark='.'; end
    if prm.showf(4)
        plot(time,res(:,n),'.','color',[0.8,0.8,0.8],'marker',mark,'markersize',prm.psize);
    end
    if prm.ptype>1
        plot(time(j),res(j,n),'b.','marker',mark,'markersize',prm.psize)
    end
    plot(time(k),res(k,n),'m.','marker',mark,'markersize',prm.psize)
    if prm.ptype==1|prm.ptype==3
        plotsep(time(j),time(j),res(j,n),'b-',prm.tint);
    end
    if prm.showf(1)
        s=sprintf('MEAN : %.4fm RMS : %.4fm',meann(res(j,n)),rmsn(res(j,n)));
        ggt('mtext',['res_',num2str(n)],s,1);
    end
    if any(prm.showf(2:3))
        t=prm.time(1)/3600:prm.time(end)/3600+1; rsm=repmat(nan,length(t),1); rss=rsm;
        for m=1:length(t)
            i=find(t(m)-0.5<=time(j)&time(j)<t(m)+0.5);
            rsm(m)=meann(res(j(i),n));
            rss(m)=stdn(res(j(i),n));
        end
        if prm.showf(2), plot(t,rsm,'r-'), plot(t,rsm,'r.'), end
        if prm.showf(3), plot(t,rss,'m-'), plot(t,rss,'m.'), end
    end
    ylabel(label{n})
    if n==1
        ti=sprintf('Postfit Residuals %s-%s : %s',prm.sat,rcv,pstr(prm.td,prm.time));
        title(ti);
    end
end

% plot residuals by Az/El ------------------------------------------------------
function PlotResidualAzEl
data=get(gcf,'userdata'); prm=data.prm;
if isempty(data.t), return, end
if strcmp(prm.rcv,'ALL'), rcv=prm.rcvs{1}; else rcv=prm.rcv; end
n=find(strcmp(rcv,prm.rcvs)); if isempty(n), return, end
i=find(any(data.index(:,3:4)==n,2));
time=data.t(i)/3600;
index=data.index(i,:);
res=data.ress(i,2);
out=data.outs(i);
azs=data.azels(i,1)*180/pi;
els=data.azels(i,2)*180/pi;
i=ones(length(res),1);
if ~strcmp(prm.sat,'ALL')
    n=find(strcmp(prm.sat,prm.sats)); if isempty(n), return, end
    i=i&any(index(:,1:2)==n,2);
end
j=find(i&~out); k=find(i&out);
clf
margin=[0.08,0.03,0.04,0.025];
ggt('subplotv','res',1,2,'move',[0,0,0,1],'link',[0,0,0,1],...
    'margin',margin,'xlim',[-180,180],'ylim',prm.range,'xtick',-180:30:180);
ylabel('Azimuth(deg)-Residuals(m)');
if prm.ptype==4, mark='o'; else mark='.'; end
if prm.showf(4)
    plot(azs,res,'.','color',[0.8,0.8,0.8],'marker',mark,'markersize',prm.psize);
end
plot(azs(j),res(j),'b.','marker',mark,'markersize',prm.psize);
plot(azs(k),res(k),'m.','marker',mark,'markersize',prm.psize);
if prm.showf(1)
    s=sprintf('MEAN : %.4fm RMS : %.4fm',meann(res(j)),rmsn(res(j)));
    ggt('mtext','res_1',s,1);
end
az=-180:prm.mesh(1):180; rsm=repmat(nan,length(az),1); rss=rsm;
for n=1:length(az)
    i=find(az(n)-prm.mesh(1)/2<=azs(j)&azs(j)<az(n)+prm.mesh(1)/2);
    rsm(n)=meann(res(j(i)));
    rss(n)=stdn(res(j(i)));
end
if prm.showf(2), plot(az,rsm,'r-'), plot(az,rsm,'r.'); end
if prm.showf(3), plot(az,rss,'m-'), plot(az,rss,'m.'); end
ti=sprintf('Postfit Residuals %s-%s : %s',prm.sat,rcv,pstr(prm.td,prm.time));
title(ti);

ggt('subplotv','res',2,2,'move',[0,0,0,1],'link',[0,0,0,1],...
    'margin',margin,'xlim',[0,90],'ylim',prm.range,'xtick',0:10:90);
ylabel('Elevation(deg)-Residuals(m)'); 
if prm.showf(3)
    plot(els,res,'.','color',[0.7,0.7,0.7],'marker',mark,'markersize',prm.psize);
end
plot(els(j),res(j),'b.','marker',mark,'markersize',prm.psize)
plot(els(k),res(k),'m.','marker',mark,'markersize',prm.psize)
if prm.showf(1)
    s=sprintf('MEAN : %.4fm RMS : %.4fm',meann(res(j)),rmsn(res(j)));
    ggt('mtext','res_2',s,1);
end
el=0:prm.mesh(2):90; rsm=repmat(nan,length(el),1); rss=rsm;
for n=1:length(el)
    i=find(el(n)-prm.mesh(2)/2<=els(j)&els(j)<el(n)+prm.mesh(2)/2);
    rsm(n)=meann(res(j(i)));
    rss(n)=stdn(res(j(i)));
end
if prm.showf(2), plot(el,rsm,'r-'), plot(el,rsm,'r.'), end
if prm.showf(3), plot(el,rss,'m-'), plot(el,rss,'m.'), end

% plot residuals by Heading/Nadir ----------------------------------------------
function PlotResidualThPh
data=get(gcf,'userdata'); prm=data.prm;
if isempty(data.t), return, end
if strcmp(prm.sat,'ALL'), sat=prm.sats{1}; else sat=prm.sat; end
n=find(strcmp(sat,prm.sats)); if isempty(n), return, end
i=find(any(data.index(:,1:2)==n,2));
time=data.t(i)/3600;
index=data.index(i,:);
res=data.ress(i,2);
out=data.outs(i);
ths=data.thphs(i,1)*180/pi;
phs=data.thphs(i,2)*180/pi;
i=ones(length(res),1);
if ~strcmp(prm.rcv,'ALL')
    n=find(strcmp(prm.rcv,prm.rcvs)); if isempty(n), return, end
    i=i&any(index(:,3:4)==n,2);
end
j=find(i&~out); k=find(i&out);
clf
margin=[0.08,0.03,0.04,0.025];
if prm.ptype==4, mark='o'; else mark='.'; end
ggt('subplotv','res',1,2,'move',[0,0,0,1],'link',[0,0,0,1],...
    'margin',margin,'xlim',[-180,180],'ylim',prm.range,'xtick',-180:30:180);
ylabel('Heading(deg)-Residuals(m)');
if prm.showf(4)
    plot(phs,res,'.','color',[0.8,0.8,0.8],'marker',mark,'markersize',prm.psize);
end
plot(phs(j),res(j),'b.','marker',mark,'markersize',prm.psize);
plot(phs(k),res(k),'m.','marker',mark,'markersize',prm.psize);
if prm.showf(1)
    s=sprintf('MEAN : %.4fm RMS : %.4fm',meann(res(j)),rmsn(res(j)));
    ggt('mtext','res_1',s,1);
end
ph=-180:prm.mesh(1):180; rsm=repmat(nan,length(ph),1); rss=rsm;
for n=1:length(ph)
    i=find(ph(n)-prm.mesh(1)/2<=phs(j)&phs(j)<ph(n)+prm.mesh(1)/2);
    rsm(n)=meann(res(j(i)));
    rss(n)=stdn(res(j(i)));
end
if prm.showf(2), plot(ph,rsm,'r-'), plot(ph,rsm,'r.'); end
if prm.showf(3), plot(ph,rss,'m-'), plot(ph,rss,'m.'); end
ti=sprintf('Postfit Residuals %s-%s : %s',sat,prm.rcv,pstr(prm.td,prm.time));
title(ti);

ggt('subplotv','res',2,2,'move',[0,0,0,1],'link',[0,0,0,1],...
    'margin',margin,'xlim',[0,15],'ylim',prm.range,'xtick',0:1:15);
ylabel('Nadir(deg)-Residuals(m)'); 
if prm.showf(4)
    plot(ths,res,'.','color',[0.7,0.7,0.7],'marker',mark,'markersize',prm.psize);
end
plot(ths(j),res(j),'b.','marker',mark,'markersize',prm.psize)
plot(ths(k),res(k),'m.','marker',mark,'markersize',prm.psize)
if prm.showf(1)
    s=sprintf('MEAN : %.4fm RMS : %.4fm',meann(res(j)),rmsn(res(j)));
    ggt('mtext','res_2',s,1);
end
th=0:prm.mesh(2)/5:15; rsm=repmat(nan,length(th),1); rss=rsm;
for n=1:length(th)
    i=find(th(n)-prm.mesh(2)/10<=ths(j)&ths(j)<th(n)+prm.mesh(2)/10);
    rsm(n)=meann(res(j(i)));
    rss(n)=stdn(res(j(i)));
end
if prm.showf(2), plot(th,rsm,'r-'), plot(th,rsm,'r.'), end
if prm.showf(3), plot(th,rss,'m-'), plot(th,rss,'m.'), end

% plot residuals in sky map ---------------------------------------------------
function PlotResidualSky
data=get(gcf,'userdata'); prm=data.prm;
if isempty(data.t), return, end
if strcmp(prm.rcv,'ALL'), rcv=prm.rcvs{1}; else rcv=prm.rcv; end
n=find(strcmp(rcv,prm.rcvs)); if isempty(n), return, end
i=find(any(data.index(:,3:4)==n,2));
index=data.index(i,:);
res=data.ress(i,2);
out=data.outs(i);
azels=data.azels(i,:);
i=ones(length(res),1);
if ~strcmp(prm.sat,'ALL')
    n=find(strcmp(prm.sat,prm.sats)); if isempty(n), return, end
    i=i&any(index(:,1:2)==n,2);
end
j=find(i&~out);
clf, hold on
[fname,fsize]=gut('getfont','g');
if prm.ptype==4, mark='o'; else mark='.'; end
set(gca,'position',[0.06,0.05,0.88,0.88],'fontname',fname,'fontsize',fsize)
ggt('skymap','fontname',fname,'fontsize',fsize);
cs=get(gcf,'colormap');
s=size(cs,1)/(prm.range(2)-prm.range(1));
if prm.showf(4)
    ggt('skyplot',data.azels,mark,'markersize',prm.psize,'color',[0.8,0.8,0.8]);
end
for n=j'
    m=floor(s*(res(n)-prm.range(1)));
    if m<1, m=1; elseif m>size(cs,1), m=size(cs,1); end
    ggt('skyplot',azels(n,:),mark,'markersize',prm.psize,'color',cs(m,:));
end
ggt('colorbarv',[0.85,0.05,0.02,0.88],prm.range,'(m)','fontname',fname,'fontsize',fsize)
ti=sprintf('Postfit Residuals %s-%s : %s',prm.sat,rcv,pstr(prm.td,prm.time));
title(ti);

% plot residuals in 3d sky map -------------------------------------------------
function PlotResidualSky3D
data=get(gcf,'userdata'); prm=data.prm;
if isempty(data.t), return, end
if strcmp(prm.rcv,'ALL'), rcv=prm.rcvs{1}; else rcv=prm.rcv; end
n=find(strcmp(rcv,prm.rcvs)); if isempty(n), return, end
i=find(any(data.index(:,3:4)==n,2));
index=data.index(i,:); time=data.t(i); res=data.ress(i,2); azels=data.azels(i,:);
clf, hold on
[fn,fs]=gut('getfont','g');
set(gca,'position',[0.06,0.05,0.88,0.88],'fontname',fn,'fontsize',fs)
r=max(abs(prm.range)); range=[-r,0,r]; scale=20/sin(prm.vangle(2)*pi/180);
ggt('skymap','dirs','rotate',prm.vangle(1),'fontname',fn,'fontsize',fs);
set(gca,'dataaspectratio',[1,1/sin(prm.vangle(2)*pi/180),1])
if prm.ptype==4, mark='o'; else mark='.'; end
if strcmp(prm.sat,'ALL'), sats=1:length(prm.sats);
else sats=find(strcmp(prm.sat,prm.sats)); end
for n=sats
    i=find(any(index(:,1:2)==n,2));
    j=find(diff(time(i))>=prm.tint*2);
    for k=[1,j'+1;j',length(i)]
        p=i(k(1):k(2));
        x=(90-azels(p,2)*180/pi).*sin(azels(p,1)+prm.vangle(1)*pi/180);
        y=(90-azels(p,2)*180/pi).*cos(azels(p,1)+prm.vangle(1)*pi/180);
        z=res(p)*scale/(range(2)-range(1));
        if any(prm.ptype==[1,3])&prm.ptype~=4
            plot(x,y+z,'-','clipping','off');
        end
        if any(prm.ptype==2:4)
            plot(x,y+z,mark,'markersize',prm.psize,'clipping','off');
        end
        plot(x,y,':','color','b');
    end
end
x=80; y=55+scale*(range-range(1))/(range(end)-range(1));
plot([x,x],[y(1),y(end)],'k','clipping','off');
for n=1:length(y)
    plot([x-2,x+2],[y(n),y(n)],'k','clipping','off');
    s=['  ',num2str(range(n))]; if range(n)==0, s=[s,' (m)']; end
    text(x+2,y(n),s,'fontname',fn,'fontsize',fs)
end
if prm.showf(1)
    s=sprintf(' MEAN : %.4fm\n RMS : %.4fm',meann(res),rmsn(res));
    text(-90,85,s,'fontname',fn,'fontsize',fs);
end
ti=sprintf('Postfit Residuals %s-%s : %s',prm.sat,rcv,pstr(prm.td,prm.time));
title(ti);

% plot residuals in sky map average/std dev ------------------------------------
function PlotResidualSkyAve(opt)
data=get(gcf,'userdata'); prm=data.prm;
if isempty(data.t), return, end
if strcmp(prm.rcv,'ALL'), rcv=prm.rcvs{1}; else rcv=prm.rcv; end
n=find(strcmp(rcv,prm.rcvs)); if isempty(n), return, end
i=find(any(data.index(:,3:4)==n,2));
index=data.index(i,:);
res=data.ress(i,2);
out=data.outs(i);
azels=data.azels(i,:)*180/pi;
i=ones(length(res),1);
if ~strcmp(prm.sat,'ALL')
    n=find(strcmp(prm.sat,prm.sats)); if isempty(n), return, end
    i=i&any(index(:,1:2)==n,2);
end
j=find(i&~out);
clf, hold on
[fname,fsize]=gut('getfont','g');
set(gca,'position',[0.06,0.05,0.88,0.88],'fontname',fname,'fontsize',fsize)
cs=get(gcf,'colormap');
s=size(cs,1)/(prm.range(2)-prm.range(1));
da=prm.mesh(1); de=prm.mesh(2);
for az=-180+da/2:da:180-da/2
    for el=de/2:de:90-de/2
        i=find(az-da/2<=azels(j,1)&azels(j,1)<az+da/2&...
               el-de/2<=azels(j,2)&azels(j,2)<el+de/2);
        if ~isempty(i)
            if opt, resm=stdn(res(j(i))); else resm=meann(res(j(i))); end
            m=floor(s*(resm-prm.range(1)));
            if m<1, m=1; elseif m>size(cs,1), m=size(cs,1); end
            ggt('skycpatch',[az,el]*pi/180,da*pi/180,de*pi/180,cs(m,:));
        else
            ggt('skycpatch',[az,el]*pi/180,da*pi/180,de*pi/180,'w');
        end
    end
end
ggt('skymap','noback','fontname',fname,'fontsize',fsize);
ggt('colorbarv',[0.85,0.05,0.02,0.88],prm.range,'(m)','fontname',fname,'fontsize',fsize)
if opt, s='Std Dev '; else s='Average '; end
ti=sprintf('Postfit Residuals %s %s-%s : %s',s,prm.sat,rcv,pstr(prm.td,prm.time));
title(ti);

% summary statistics -----------------------------------------------------------
function SummaryStats
data=get(gcf,'userdata'); prm=data.prm;
if isempty(data.t), return, end
s={};
s=adds(s,'                        ESTIMATION STATISTICS');
s=adds(s,'TIME    : %s GPST',pstr(prm.td,prm.time));
s=adds(s,'');
s=AddStats(s,prm.sats,data.t,data.index,data.ress,data.outs,0);
s=adds(s,'');
s=AddStats(s,prm.rcvs,data.t,data.index,data.ress,data.outs,1);
if isempty(prm.files), ti=prm.dirs.est; else ti=prm.files{1}; end
Viewer(s,ti);

% sat antenna pcv ---------------------------------------------------------------
function SatAntPcv
data=get(gcf,'userdata'); prm=data.prm;
if isempty(data.t), return, end
s={};
s=adds(s,'function prm=prm_satpcv');
s=adds(s,'prm={');
for n=1:length(prm.sats)
    i=find(any(data.index(:,1:2)==n,2));
    time=data.t(i)/3600;
    res=data.ress(i,2);
    out=data.outs(i);
    ths=data.thphs(i,1)*180/pi;
    th=0:15; pcv=zeros(1,length(th));
    for m=1:length(th)
        i=find(th(m)-0.5<=ths&ths<th(m)+0.5&~out);
        if ~isempty(i), pcv(m)=meann(res(i)); end
    end
    strs=''; for m=1:length(th), strs=[strs,sprintf(',%5.1f',pcv(m)*1000)]; end
    s=adds(s,['''',prm.sats{n},'''',strs]);
end
s=adds(s,'};');
s=adds(s,'for n=1:size(prm,1), for m=2:17, prm{n,m}=prm{n,m}*1E-3; end, end');
Viewer(s,'Satellite Antenna PCV Parameters');

% display estimation statistics ------------------------------------------------
function s=AddStats(s,name,t,index,ress,outs,rs)
if rs, type='receiver'; else type=' sat   '; end
s=adds(s,'                                             residuals (m)');
s=adds(s,'%s   total obs  outlier prefit-rms postfit-rms ph.b-rms   ph.b-ave ',type);
s=adds(s,'-------------------------------------------------------------------------');
nobs=zeros(length(name),1); nout=nobs; rmsp=nobs; rmsf=nobs; rmsb=nobs; aveb=nobs;
for n=1:length(name)
    if rs, j=find(index(:,3)==n|index(:,4)==n);
    else   j=find(index(:,1)==n|index(:,2)==n); end
    if ~isempty(j)
        nobs(n)=length(j);
        k=j(find(~outs(j)));
        nout(n)=nobs(n)-length(k);
        if ~isempty(k)
            rmsp(n)=sum(ress(k,1).^2);
            rmsf(n)=sum(ress(k,2).^2);
            rmsb(n)=sum(ress(k,3).^2);
            aveb(n)=sum(ress(k,3));
        end
    end
    s=adds(s,'%-7s : %9d %6d %10.4f %10.4f %10.4f %10.4f',name{n},nobs(n),nout(n),...
           rms(rmsp(n),nobs(n)),rms(rmsf(n),nobs(n)),rms(rmsb(n),nobs(n)),...
           ave(aveb(n),nobs(n)));
end
s=adds(s,'-------------------------------------------------------------------------');
s=adds(s,'TOTAL   : %9d %6d %10.4f %10.4f %10.4f %10.4f',sum(nobs),sum(nout),...
       rms(sum(rmsp),sum(nobs)),rms(sum(rmsf),sum(nobs)),rms(sum(rmsb),sum(nobs)),...
       ave(sum(aveb),sum(nobs)));

% calculate satellite/receiver directions/eclipse -------------------------------
function [azels,thphs,ecl]=SatDir(td,time,index,sats,rcvs,dirs,srcnav)
azels=[]; thphs=[]; ecl=[];
if isempty(time)|isempty(index), return, end
[nav,inav]=readnav(td,time,sats,rcvs,dirs.nav,srcnav);
posr=shiftdim(readpos(td,time(1),rcvs,'','approx'),1);
[t,it,jt]=unique(time);
poss=zeros(3,length(t),length(sats)); possun=zeros(3,length(t));
for n=1:length(sats)
    navs=nav(find(inav==n),:);
    for m=1:length(t)
        poss(:,m,n)=navtostate(td,t(m),navs);
    end
end
for n=1:length(t)
    tut=td+(t(n)-19)/86400; U=ecsftoecef(tut);
    [rsun,rmoon]=sunmoonpos(tut); possun(:,n)=U*rsun; posmoon(:,n)=U*rmoon;
end
azels=zeros(length(time),2);
thphs=zeros(length(time),2);
for n=1:length(time)
    i=index(n,1); j=index(n,3); k=jt(n);
    azels(n,:)=satazel(poss(:,k,i),posr(:,j));
    thphs(n,:)=rcvdir(poss(:,k,i),posr(:,j),possun(:,k));
    ecl(n,1)=~isnan(poss(1,k,i))&shadowfunc(poss(:,k,i),possun(:,k),posmoon(:,k))==0;
end

% receiver directions in satellite fixed coordinate ----------------------------
function thph=rcvdir(psat,prcv,psun)
[ex,ey,ez]=satfixed(psat,psun);
er=[ex,ey,ez]'*(prcv-psat); er=er/norm(er);
thph=[acos(er(3)),atan2(er(2),er(1))]; % [theta,phi]

% read data --------------------------------------------------------------------
function ReadData
h=gcf;
data=get(h,'userdata'); prm=data.prm;
[prm.td,ts]=caltomjd(prm.tstart); [tn,te]=caltomjd(prm.tend);
prm.time=ts:prm.tint:te+(tn-prm.td)*86400; if isempty(prm.time), return, end
set(h,'pointer','watch');
[t,index,ress,outs,prm.sats]=readstats(prm.td,prm.time,prm.dirs.est,prm.fb,...
                                       prm.tunit,prm.rcvs);
if ~isempty(ress)
    i=find(~isnan(ress(:,2)));
    data.t=t(i); data.index=index(i,:); data.ress=ress(i,:); data.outs=outs(i);
    [data.azels,data.thphs,data.ecl]=SatDir(prm.td,data.t,data.index,prm.sats,...
                                            prm.rcvs,prm.dirs,prm.srcnav);
end
set(h','name',['Residuals/Statistics : ',prm.dirs.est,' (',prm.fb,')']);
data.prm=prm; set(h,'userdata',data);
set(h,'pointer','arrow');

delete(gut('geth','sats'));
delete(gut('geth','rcvs'));
gut('newmenu','sats','&Satellite',{'ALL',prm.sats{:}},{'ALL',prm.sats{:}},[mfilename,' cbSat']);
gut('newmenu','rcvs','&Receiver',{'ALL',prm.rcvs{:}},{'ALL',prm.rcvs{:}},[mfilename,' cbRcv']);
gut('chkmenu','sats','ALL');
gut('chkmenu','rcvs','ALL');

% plot separated ---------------------------------------------------------------
function plotsep(t,x,y,opt,tint)
i=1; for j=find(diff(t)>=tint*2/3600)', plot(x(i:j),y(i:j),opt); i=j+1; end
plot(x(i:end),y(i:end),opt);

% average/rms ------------------------------------------------------------------
function r=ave(sums,no), if no>0, r=sums/no; else r=0; end
function r=rms(sums,no), if no>0, r=sqrt(sums/no); else r=0; end

% mean/std without nan ---------------------------------------------------------
function m=meann(x)
for n=1:size(x,2)
    i=find(~isnan(x(:,n)));
    if isempty(i), m(1,n)=nan; else m(1,n)=mean(x(i,n),1); end
end

function s=stdn(x)
for n=1:size(x,2)
    i=find(~isnan(x(:,n)));
    if isempty(i), s(1,n)=nan; else s(1,n)=std(x(i,n),1); end
end

function s=rmsn(x)
for n=1:size(x,2)
    i=find(~isnan(x(:,n)));
    if isempty(i), s(1,n)=nan; else s(1,n)=sqrt(mean(x(i,n).^2,1)); end
end

% show viewer ------------------------------------------------------------------
function Viewer(str,ti)
h=gut('newviewer','',['Summary of Statistics : ',ti],[600,419,0,-48],'str',str);

% date/time string -------------------------------------------------------------
function s=tstr(td,ts)
s=sprintf('%04d/%02d/%02d %02d:%02d:%02.0f',mjdtocal(td,ts));

function s=pstr(td,time)
s=[tstr(td,time(1)),'-',tstr(td,time(end))];

% add string -------------------------------------------------------------------
function s=adds(s,varargin), s={s{:},sprintf(varargin{:})};
