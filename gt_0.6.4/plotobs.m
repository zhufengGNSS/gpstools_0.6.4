function h=plotobs(varargin)
%-------------------------------------------------------------------------------
% [system] : GpsTools
% [module] : plot observation data
% [func]   : plot observation data
% [argin]  : 'prm',prm   = parameters
% [argout] : h = figure handle
% [note]   :
% [version]: $Revision: 12 $ $Date: 2008-11-25 10:02:15 +0900 (火, 25 11 2008) $
%            Copyright(c) 2004-2006 by T.Takasu, all rights reserved
% [history]: 04/11/18  0.1  new
%            04/12/16  0.2  support gui interface
%-------------------------------------------------------------------------------
if nargin>0&strncmp(varargin{1},'cb',2), feval(varargin{:});
else h=PlotMain(varargin{:}); end

% main -------------------------------------------------------------------------
function h=PlotMain(varargin)
h=gut('newfigg','','Observation Data',[600,400,0,-68],[mfilename,' cbClose']);
if isempty(h), return, end
prm=loadprm('prm_plotobs','prm_plotobs_def');
if nargin>=2&strcmp(varargin{1},'prm'), prm=varargin{2}; end
menu1={'Data Availability','Reset Time Range','-',...
       'L1','L2','P1','P2','L1-P1','L2-P2','L1-TD','L2-TD','-','LC','LG','LG-D',...
       'PC','PG','S1','S2','SC','SG','WL','NL','MW','WL-TD','NL-TD','LG+PG','-',...
       'MP1','MP2','MP1+MP2','MP1-MP2','-','Timetag Offset','Clock Bias','STEC','VTEC','Azimuth','Elevation'};
gut('newmenu','data','&Data ',...
    {'&Read...','&New Window','-','&Export Plot...','-','&Raw Obs Data...',...
     '&TEQC...','-','&Options...'},{},...
    {[mfilename,' cbRead'],[mfilename,' cbNew'],'',[mfilename,' cbExport'],'',...
     [mfilename,' cbView'],[mfilename,' cbTeqc'],'',[mfilename,' cbOpt']});
gut('newmenu','plot','&Plot',menu1,menu1,[mfilename,' cbPlot']);
gut('newmenu','sats','&Satellite',prm.sats,prm.sats,[mfilename,' cbSat']);
gut('newmenu','rcvs','&Receiver',prm.rcvs,prm.rcvs,[mfilename,' cbRcv']);
gut('chkmenu','plot',prm.type);
gut('chkmenu','sats',prm.sat);
gut('chkmenu','rcvs',prm.rcv);
set(gut('newbtnh','sbtn',[0,0,102,16],{'','','<','>','<','>'},[mfilename,' cbSwBtn']),...
    'fontsize',8);
set(gcf,'resizefcn',[mfilename,' cbResize']);
data.tspan=[0,24]; data.zo=[]; data.izo=[]; data.zc=[]; data.izc=[];
prm.td=caltomjd([2000,1,1]); prm.time=0; prm.tspan=[0,1];
data.rpos=[]; data.anttype={''}; data.rcvtype={''}; data.rstat=[]; data.azel=[];
data.slip=[]; data.prm=prm;
set(gcf,'userdata',data)
LocatePos;
UpdateBtns;
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
sel1={{'RINEX','RINEX(3H)','RINEX(1H)','RINEX(15min)'},...
      {'rinex','rinex3','rinex1','rinexh'}};
sel2={{'Combined','RINEX','RINEX(3H)','RINEX(1H)','RINEX(15min)'},...
      {'brdc','rinex','rinex3','rinex1','rinexh'}};
data=get(gcf,'userdata'); prm=data.prm;
prm1={
't','Start Time (GPST)',[]
't','End Time (GPST)',  []
'n','Observation Unit Time (hr)',24
'b','Satellites',[mfilename,' cbSelSat']
'b','Receiver',[mfilename,' cbSelRcv']
'p','Raw Observation Data', sel1
'p','Navigation Messages',  sel2
};
prm2={
' ','Observation Data Directory',''
'd','',prm.dirs.obs
' ','Clean Observation Data Directory',''
'd','',prm.dirs.obc
' ','Navigation Message Directory',''
'd','',prm.dirs.nav
};
q='';
gut('newdlg','','Read Data',[314,331]);
gut('newprms','prm1',[12,150,295,25,118],prm1);
gut('newprms','prm2',[12,34,295,18,118],prm2);
gut('newokcancelbtn','',[127,4,180,23]);
gut('setprms','prm1',prm.tstart,prm.tend,prm.tunit,'','',prm.src,prm.navsrc);
gut('setudata','prm1_4',prm.sats);
gut('setudata','prm1_5',prm.rcvs);
if ~gut('waitok'), return, end
[prm.tstart,prm.tend,prm.tunit,q,q,prm.src,prm.navsrc]=gut('getprms','prm1');
[q,prm.dirs.obs,q,prm.dirs.obc,q,prm.dirs.nav]=gut('getprms','prm2');
prm.sats=gut('getudata','prm1_4');
prm.rcvs=gut('getudata','prm1_5');
if ~isempty(prm.rcvs), prm.rcv=prm.rcvs{1}; end
if ~isempty(prm.sats), prm.sat=prm.sats{1}; end
close
prm.files={};
delete(gut('geth','sats'));
delete(gut('geth','rcvs'));
gut('newmenu','sats','&Satellite',prm.sats,prm.sats,[mfilename,' cbSat']);
gut('newmenu','rcvs','&Receiver',prm.rcvs,prm.rcvs,[mfilename,' cbRcv']);
gut('chkmenu','sats',prm.sat);
gut('chkmenu','rcvs',prm.rcv);
data.prm=prm; set(gcf,'userdata',data)
SaveSetting; ReadData; UpdatePlot;

% callback on menu select receiver --------------------------------------------
function cbSelRcv
[rcvs,ok]=editlist('Receiver',get(gcbo,'userdata'),'rcv');
if ok, set(gcbo,'userdata',rcvs); end

% callback on menu select satellite -------------------------------------------
function cbSelSat
[sats,ok]=editlist('Satellites',get(gcbo,'userdata'));
if ok, set(gcbo,'userdata',sats); end

% callback on menu export plot -------------------------------------------------
function cbExport
persistent file dirs
if isempty(file), file='plotobs'; end, if isempty(dirs), dirs=''; end
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
if isempty(strmatch('Data Availability',prm.type))
    sat=prm.sat;
    for n=1:length(prm.sats)
        prm.sat=prm.sats{n}; data.prm=prm; set(f,'userdata',data)
        figure(f), UpdatePlot;
        print(format{:},opts{:},sprintf('%s%02d',fs,n));
    end
    prm.sat=sat; data.prm=prm; set(f,'userdata',data); UpdatePlot;
else
    rcv=prm.rcv;
    for n=1:length(prm.rcvs)
        prm.rcv=prm.rcvs{n}; data.prm=prm; set(f,'userdata',data)
        figure(f), ReadData; UpdatePlot;
        print(format{:},opts{:},sprintf('%s%02d',fs,n));
    end
    prm.rcv=rcv; data.prm=prm; set(f,'userdata',data); ReadData; UpdatePlot;
end

% callback on menu view --------------------------------------------------------
function cbView, ViewObs

% callback on menu teqc --------------------------------------------------------
function cbTeqc, TeqcObs

% callback on menu options -----------------------------------------------------
function cbOpt
sel1={{'OFF','ON'},0:1};
sel2={{'Dot','Line','Line and Dot'},1:3};
sel3={{'Raw Obs Data','Clean Obs Data','Both Obs Data'},1:3};
data=get(gcf,'userdata'); prm=data.prm;
prm1={
'p','Plot Obs Data',                   sel3
'p','Eliminate Phase Bias',            sel1
'p','Show Number of Visible Sats.',    sel1
'p','Show Antenna/Receiver Type',      sel1
'p','Plot Style',                      sel2
'c','Plot Color of Raw Obs. Data',     nan
'c','Plot Color of Clean Obs. Data',   nan
'c','Plot Color of Cycle-Slip',        nan
'n','Max Arc Gap (sec)',               60
's','Value Range',                     [-inf,inf]
};
prm2={
'e','TEQC Options',                    prm.teqcopt
};
q='';
gut('newdlg','','Options',[320,315]);
gut('newprms','prm1',[12,60,300,25,120],prm1);
gut('newprms','prm2',[12,35,300,25,220],prm2);
gut('setprms','prm1',prm.plotf,prm.nobias,prm.showf(1),prm.showf(3),prm.ptype,...
    prm.color{1},prm.color{2},prm.color{3},prm.mgap,prm.range);
gut('newokcancelbtn','',[135,4,180,23]);
if ~gut('waitok'), return, end
[prm.plotf,prm.nobias,prm.showf(1),prm.showf(3),prm.ptype,prm.color{1},prm.color{2},...
 prm.color{3},prm.mgap,range]=gut('getprms','prm1');
prm.teqcopt=gut('getprms','prm2');
if range(1)<range(2), prm.range=range; else prm.range=[-inf,inf]; end
close
data.prm=prm; set(gcf,'userdata',data);
SaveSetting; UpdatePlot;

% callback on menu plot --------------------------------------------------------
function cbPlot
data=get(gcf,'userdata'); prm=data.prm;
switch get(gcbo,'label')
case 'Reset Time Range'
    data.tspan=prm.tspan;
    set(gcf,'userdata',data)
    UpdatePlot
otherwise
    gut('togglechk',gcbo);
    prm.type=gut('getchk','plot');
    data.prm=prm; set(gcf,'userdata',data)
    UpdatePlot
end

% callback on menu satellite ---------------------------------------------------
function cbSat
data=get(gcf,'userdata'); prm=data.prm;
prm.sat=get(gcbo,'userdata');
gut('unchkmenu','sats');
gut('chkmenu','sats',prm.sat);
data.prm=prm; set(gcf,'userdata',data);
UpdatePlot

% callback on menu receiver ----------------------------------------------------
function cbRcv
data=get(gcf,'userdata'); prm=data.prm;
prm.rcv=get(gcbo,'userdata');
gut('unchkmenu','rcvs');
gut('chkmenu','rcvs',prm.rcv);
data.prm=prm; set(gcf,'userdata',data);
ReadData
UpdatePlot

% callback on change axis ------------------------------------------------------
function cbAxis
data=get(gcf,'userdata');
data.tspan=get(gca,'xlim');
set(gcf,'userdata',data)

% callback on push swithch button ----------------------------------------------
function cbSwBtn
data=get(gcf,'userdata'); prm=data.prm;
i=find(strcmp(prm.sats,prm.sat));
j=find(strcmp(prm.rcvs,prm.rcv));
switch get(gcbo,'tag')
case 'sbtn_1', if prm.plotf==3, prm.plotf=1; else prm.plotf=prm.plotf+1; end
case 'sbtn_2', if prm.showf(2), prm.showf(2)=0; else prm.showf(2)=1; end
case 'sbtn_3', prm.sat=prm.sats{max(i-1,1)};
case 'sbtn_4', prm.sat=prm.sats{min(i+1,length(prm.sats))};
case 'sbtn_5', prm.rcv=prm.rcvs{max(j-1,1)};
case 'sbtn_6', prm.rcv=prm.rcvs{min(j+1,length(prm.rcvs))};
end
gut('unchkmenu','sats'); gut('chkmenu','sats',prm.sat);
gut('unchkmenu','rcvs'); gut('chkmenu','rcvs',prm.rcv);
data.prm=prm; set(gcf,'userdata',data);
if any(strcmp(get(gcbo,'tag'),{'sbtn_5','sbtn_6'})), ReadData; end
UpdateBtns;
UpdatePlot;

% locate position --------------------------------------------------------------
function LocatePos
p=get(gcf,'position'); m=15;
gut('setpos','sbtn_1',[m,p(4)-16]);
gut('setpos','sbtn_2',[m+16,p(4)-16]);
gut('setpos','sbtn_3',[p(3)-m-64,p(4)-16]);
gut('setpos','sbtn_4',[p(3)-m-48,p(4)-16]);
gut('setpos','sbtn_5',[p(3)-m-32,p(4)-16]);
gut('setpos','sbtn_6',[p(3)-m-16,p(4)-16]);

% save settingsrs --------------------------------------------------------------
function SaveSetting
data=get(gcf,'userdata'); prm=data.prm;
[path,file]=fileparts(which(mfilename));
save(fullfile(path,'settings','prm_plotobs.mat'),'prm');

% update switch button ---------------------------------------------------------
function UpdateBtns
data=get(gcf,'userdata'); prm=data.prm;
ls1={'o','c','b'}; ls2={' ','s'};
gut('setstring','sbtn_1',ls1{prm.plotf});
gut('setstring','sbtn_2',ls2{prm.showf(2)+1});

% read observation data --------------------------------------------------------
function ReadData
data=get(gcf,'userdata'); prm=data.prm;
[prm.td,ts]=caltomjd(prm.tstart); [tn,te]=caltomjd(prm.tend);
prm.time=ts:prm.tint:te+(tn-prm.td)*86400;
if ~isempty(prm.time)
    prm.tspan=[prm.time(1),prm.time(end)]/3600; data.tspan=prm.tspan;
    h=gcf; set(h,'pointer','watch');
    [data.zo,data.izo,data.rpos]=...
        readobs(prm.td,prm.time,prm.sats,prm.rcv,prm.dirs.obs,prm.src,'',inf);
    [data.zc,data.izc,rp,ad,data.anttype,data.rcvtype,data.rstat,data.azel,data.slip]=...
        readobs(prm.td,prm.time,prm.sats,prm.rcv,prm.dirs.obc,'obsc','',inf,prm.tunit);
    data.izo(:,3)=find(strcmp(prm.rcv,prm.rcvs));
    data.izc(:,3)=find(strcmp(prm.rcv,prm.rcvs));
    ti=[prm.dirs.obs,' (',prm.src,') : ',prm.dirs.obc];
    set(gcf,'name',['Observation Data : ',ti]);
    set(h,'pointer','arrow');
else
    prm.time=ts; prm.tspan=[0,24]; data.tspan=prm.tspan;
end
data.prm=prm; set(gcf,'userdata',data)

% update observation data plot -------------------------------------------------
function UpdatePlot
data=get(gcf,'userdata'); prm=data.prm;
clf
if any(strcmp(prm.type,'Data Availability'))
    PlotAvailability(data.izo,data.izc,data.slip,data.tspan,data.anttype{1},...
                     data.rcvtype{1},prm);
    return;
end
margin=[0.08,0.03,0.04,0.012];
for n=1:length(prm.type)
    if n<length(prm.type), topts='nolabel'; else topts=''; end
    ggt('subplotv','obs',n,length(prm.type),'move',[1,1,1,1],'link',[1,0,1,0],...
        'taxis',prm.td,'topts',topts,'margin',margin,...
        'xlim',data.tspan,'ylim',prm.range,'callback',[mfilename,' cbAxis'],...
        'callbackd',[mfilename,' cbAva'])
    zmin=[]; zmax=[]; unit=''; bias=[];
    if any(prm.plotf==[1,3])
        [z1,z2,unit]=PlotObs(1,data.zo,data.izo,data.rpos,[],[],prm.type{n},data.slip,data.tspan,prm);
        zmin=min([zmin,z1]); zmax=max([zmax,z2]);
    end
    if any(prm.plotf==[2,3])
        [z1,z2,unit]=PlotObs(2,data.zc,data.izc,data.rpos,data.rstat,data.azel,prm.type{n},data.slip,data.tspan,prm);
        zmin=min([zmin,z1]); zmax=max([zmax,z2]);
    end
    switch prm.type{n}
    case 'Timetag Offset', range=[-1E-2,1E-2];
    case 'Elevation', range=[0,90];
    case 'Azimuth', range=[-180,180];
    otherwise
        range=prm.range;
        if isinf(range(1))&~isempty(zmin), range(1)=zmin; end
        if isinf(range(2))&~isempty(zmax), range(2)=zmax; end
        if range(2)<=range(1), range=[-inf,inf]; end
    end
    if all(~isnan(range)), ylim(range); end
    ylabel([prm.type{n},' (',unit,')']);
    if n==1
        ti=[tstr(prm.td,prm.time(1)),'-',tstr(prm.td,prm.time(end))];
        title(['Observation Data ',prm.sat,'-',prm.rcv,' : ',ti])
    end
end

% plot observation data --------------------------------------------------------
function [zmin,zmax,unit]=PlotObs(obs,z,iz,rpos,rstat,azel,type,slip,tspan,prm)
zmin=[]; zmax=[]; unit='';
nbs={'P1','P2','PC','PG','LG-D','L1-TD','L2-TD','WL-TD','Timetag Offset','Clock Bias','Elevation'};
[t,z,iz,unit]=GenObs(type,z,iz,rpos,rstat,azel,prm); if isempty(t), return, end
if prm.showf(2)&~isempty(slip)
    n=find(strcmp(prm.sats,prm.sat));
    for tt=slip(slip(:,2)==n&slip(:,4)<9,1)'
        plot([tt,tt]/3600,[-1E9,1E9],'-','color',prm.color{3});
    end
    for tt=slip(slip(:,2)==n&slip(:,4)==9,1)'
        plot([tt,tt]/3600,[-1E9,1E9],'-','color','g');
    end
end
if obs==2
    gap=find(iz(:,4)==1)'-1; gap(gap<=0)=[];
else
    gap=find(diff(t)>prm.mgap)';
end
for i=[1,gap+1;gap,length(t)]
    j=i(1):i(2);
    if prm.nobias&~any(strcmp(type,nbs))
        k=j(tspan(1)*3600<=t(j)&t(j)<=tspan(2)*3600);
        z(j)=z(j)-meann(z(k));
    end
    if obs==2
        k=find(iz(j,4)==1); plot(t(j(k))/3600,z(j(k)),'o','color',prm.color{obs});
        k=find(iz(j,4)==2); plot(t(j(k))/3600,z(j(k)),'s','color',prm.color{obs});
    end
    if any(prm.ptype==[1,3]), plot(t(j)/3600,z(j),'.','color',prm.color{obs}); end
    if any(prm.ptype==[2,3])
        p=1;
        for q=find(iz(j,4)==2)'
            plot(t(j(p:q))/3600,z(j(p:q)),'-','color',prm.color{obs});
            p=q+1;
        end
        plot(t(j(p:end))/3600,z(j(p:end)),'-','color',prm.color{obs});
    end
    k=j(tspan(1)*3600<=t(j)&t(j)<=tspan(2)*3600);
    zmin=min([zmin;z(k)]);
    zmax=max([zmax;z(k)]);
end

% plot observation data availability -------------------------------------------
function PlotAvailability(izo,izc,slip,tspan,anttype,rcvtype,prm)
clf
margin=[0.08,0.03,0.04,0.012];
label={}; for n=1:length(prm.sats), label={prm.sats{n},label{:}}; end
ggt('subplotv','ava',1,1,'move',[1,0,1,0],'taxis',prm.td,'margin',margin,...
    'xlim',tspan,'ylim',[0.5,max(1,length(prm.sats))+0.5],...
    'ytick',1:length(prm.sats),'yticklabel',label,'ygrid','off',...
    'ticklength',[0,0],'callback',[mfilename,' cbAxis'],...
    'callbackd',[mfilename,' cbData'])
for n=1:length(prm.sats)
    y=length(prm.sats)-n+1;
    if any(prm.plotf==[1,3])&~isempty(izo), PlotBar(izo(izo(:,2)==n,1),y,prm.color{1}); end
    if any(prm.plotf==[2,3])&~isempty(izc), PlotBar(izc(izc(:,2)==n,1),y,prm.color{2}); end
    plot([-1E6,1E6],[y,y]+0.5,'k:')
    if prm.showf(2)&~isempty(slip)
        PlotSlip(slip(slip(:,2)==n&slip(:,4)<9,1),y,prm.color{3});
        PlotSlip(slip(slip(:,2)==n&slip(:,4)==9,1),y,'g');
    end
end
if prm.showf(3)
    ggt('mtext','ava_1',rcvtype,1);
    ggt('mtext','ava_1',anttype,2);
end
ti=[tstr(prm.td,prm.time(1)),'-',tstr(prm.td,prm.time(end))];
title(['Observation Data Availability ',prm.rcv,' : ',ti]);

if ~prm.showf(1), return, end
ggt('subplotv','',1,1,'move',[1,0,1,0],'link',[1,0,1,0],'taxis',prm.td,...
    'topts','nolabel','margin',margin,'xlim',tspan,'ylim',[0,16],...
    'box','off','color','none','yaxislocation','right','ticklength',[1E-3,0])
grid off;
if any(prm.plotf==[2,3]), PlotObsCount(izc,prm.time,'r'); end

% plot bar ---------------------------------------------------------------------
function PlotBar(t,y,color)
if isempty(t), return, end
tt=round(diff(t)); tint=max(min(tt),1); i=find(tt>tint)';
for j=[1,i+1;i,length(t)]
    plot([t(j(1)),t(j(2))]/3600,[y,y],'-','color',color,'linewidth',3)
end

% plot slip --------------------------------------------------------------------
function PlotSlip(t,y,color)
tt=[t';t';t']; tt(3,:)=nan; y=repmat(y-0.25,3,length(t)); y(1,:)=y(1,:)+0.5;
plot(tt(:)/3600,y(:),'-','color',color);

% plot obs count ---------------------------------------------------------------
function PlotObsCount(iz,time,color)
if size(iz,1)<2, return, end
t=round(iz(:,1)); tint=min(diff(unique(t)));
time=((floor(time(1)/tint):floor(time(end)/tint))*tint)';
[tt,i,j]=intersect(time,t);
nobs=zeros(length(time),1); nobs(i)=j-[0;j(1:end-1)];
plot(time/3600,nobs,'-','color',color);

% callback on double click on plot ---------------------------------------------
function cbData
data=get(gcf,'userdata'); prm=data.prm;
p=get(gca,'currentpoint'); yl=get(gca,'ylim');
n=round(length(prm.sats)*(yl(2)-p(1,2))/(yl(2)-yl(1))+0.5);
if n<1|length(prm.sats)<n, return, end
prm.sat=prm.sats{n};
prm.type(strcmp(prm.type,'Data Availability'))=[];
gut('unchkmenu','plot','Data Availability');
data.prm=prm; set(gcf,'userdata',data)
UpdatePlot

function cbAva
data=get(gcf,'userdata'); prm=data.prm;
prm.type={'Data Availability',prm.type{:}};
gut('chkmenu','plot',prm.type);
data.prm=prm; set(gcf,'userdata',data)
UpdatePlot

% generate observation combination ---------------------------------------------
function [t,zz,iz,unit]=GenObs(type,z,iz,rpos,rstat,azel,prm)
C=299792458; f1=1.57542E9; f2=1.22760E9;
lam1=C/f1; lam2=C/f2; lam4=1/(1/lam1-1/lam2); lam5=1/(1/lam1+1/lam2);
c1=f1^2/(f1^2-f2^2); c2=f2^2/(f1^2-f2^2);
t=[]; zz=[]; unit='';
n=find(strcmp(prm.sat,prm.sats));
m=find(strcmp(prm.rcv,prm.rcvs));
if isempty(iz)|isempty(n)|isempty(m), return, end
i=find(iz(:,2)==n&iz(:,3)==m); if isempty(i), return, end
tt=iz(:,1); z=z(i,:); iz=iz(i,:); t=round(iz(:,1)); zz=nan*t;
gap=find(diff(t)>prm.mgap)';
unit='m';

switch upper(type)
case 'P1', zz=z(:,3);
case 'P2', zz=z(:,4);
case 'L1', zz=lam1*z(:,1);
case 'L2', zz=lam2*z(:,2);
case 'LC', zz=c1*lam1*z(:,1)-c2*lam2*z(:,2);
case 'LG', zz=lam1*z(:,1)-lam2*z(:,2);
case 'LG-D'
    for i=[1,gap+1;gap,size(iz,1)]
        j=i(1):i(2);
        if length(j)>=2, zz(j)=[nan;diff(lam1*z(j,1)-lam2*z(j,2))./diff(iz(j,1))]; end
        unit='m/sec';
    end
case 'L1-P1', zz=lam1*z(:,1)-z(:,3);
case 'L2-P2', zz=lam2*z(:,2)-z(:,4);
case 'L1-TD'
    for i=[1,gap+1;gap,size(iz,1)]
        j=i(1):i(2);
        if length(j)>=4, zz(j)=[nan;nan;diff(diff(diff(z(j,1))./diff(iz(j,1))));nan]; end
        unit='cycle/sec';
    end
case 'L2-TD'
    for i=[1,gap+1;gap,size(iz,1)]
        j=i(1):i(2);
        if length(j)>=4, zz(j)=[nan;nan;diff(diff(diff(z(j,2))./diff(iz(j,1))));nan]; end
        unit='cycle/sec';
    end
case 'WL-TD'
    for i=[1,gap+1;gap,size(iz,1)]
        j=i(1):i(2);
        if length(j)>=4, zz(j)=[nan;nan;diff(diff(diff(z(j,1)-z(j,2))./diff(iz(j,1))));nan]; end
        unit='cycle/sec';
    end
case 'NL-TD'
    for i=[1,gap+1;gap,size(iz,1)]
        j=i(1):i(2);
        if length(j)>=4, zz(j)=[nan;nan;diff(diff(diff(z(j,1)+z(j,2))./diff(iz(j,1))));nan]; end
        unit='cycle/sec';
    end
case 'LG+PG'
    for i=[1,gap+1;gap,size(iz,1)]
        j=i(1):i(2);
        zz(j)=lam1*z(j,1)-lam2*z(j,2)+z(j,3)-z(j,4);
    end
case 'PC', zz=z(:,3:4)*[c1;-c2];
case 'PG', zz=z(:,3)-z(:,4);
case 'S1', if size(z,2)>=6, zz=z(:,6); end
case 'S2', if size(z,2)>=7, zz=z(:,7); end
case 'SC', if size(z,2)>=5, zz=z(:,5); end
case 'SG', if size(z,2)>=7, zz=z(:,6)-z(:,7); end
case 'WL', zz=lam4*(z(:,1)-z(:,2));
case 'NL', zz=lam5*(z(:,1)+z(:,2));
case 'MW', zz=lam4*(z(:,1)-z(:,2))-lam5*(z(:,3)/lam1+z(:,4)/lam2);
case {'MP1','MP2','MP1+MP2','MP1-MP2'}
    mp1=z(:,3)-(2*c2+1)*lam1*z(:,1)+2*c2*lam2*z(:,2);
    mp2=z(:,4)-2*c1*lam1*z(:,1)+(2*c1-1)*lam2*z(:,2);
    switch type
    case 'MP1', zz=mp1-meann(mp1);
    case 'MP2', zz=mp2-meann(mp2);
    case 'MP1+MP2', zz=mp1+mp2-meann(mp1+mp2);
    case 'MP1-MP2', zz=mp1-mp2-meann(mp1-mp2);
    end
case 'TIMETAG OFFSET'
    t=tt; zz=tt(:)-round(tt(:)); iz=zeros(length(t),4); unit='sec';
case 'CLOCK BIAS'
    if ~isempty(rstat)
        [tt,j,k]=intersect(t,round(rstat{1}(:,1)));
        zz(j)=rstat{1}(k,end);
    end
case 'STEC'
    lg=-(lam1*z(:,1)-lam2*z(:,2));
    zz=(lg-mean(lg-(z(:,3)-z(:,4))))*f1^2*f2^2/(f2^2-f1^2)/40.3E16;
    unit='TECU';
case 'VTEC'
    Re=6378137; H=350000;
    if ~isempty(azel)
        lg=-(lam1*z(:,1)-lam2*z(:,2));
        stec=lg-mean(lg-(z(:,3)-z(:,4)));
        cosz=sqrt(1-(Re*sin(pi/2-azel(i,2))/(Re+H)).^2);
        zz=stec.*cosz*f1^2*f2^2/(f2^2-f1^2)/40.3E16;
    end
    unit='TECU';
case 'AZIMUTH',   if ~isempty(azel), zz=azel(i,1)*180/pi; end, unit='deg';
case 'ELEVATION', if ~isempty(azel), zz=azel(i,2)*180/pi; end, unit='deg';
otherwise, disp(['waring : unknown observable : ',type])
end

% view observation file --------------------------------------------------------
function ViewObs
h=gcf; set(h,'pointer','watch');
gut('newviewer','',['Observation Data'],[600,419,0,-20],'file',GetObsFile);
set(h,'pointer','arrow');

% teqc observation file --------------------------------------------------------
function TeqcObs
data=get(gcf,'userdata'); prm=data.prm;
cmd=which('teqc.exe'); if isempty(cmd), return, end
[obs,nav]=GetObsFile; if isempty(obs), return, end
h=gcf; set(h,'pointer','watch');
obss=''; navs='';
for n=1:length(obs)
    [obsu{n},stat,msg]=uncompact(obs{n},'org','dir',pwd);
    if stat==0, obss=[obss,' "',obsu{n},'"']; end
end
for n=1:length(nav), navs=[navs,' -nav "',nav{n},'"']; end
[s,err]=dos(['"',cmd,'" ',prm.teqcopt,navs,obss]);
if s==0
    rep=obsu{end}; rep(end)='s';
    if exist(rep)
        str=textread(rep,'%s','delimiter','\n','whitespace',''); delete(rep);
        gut('newviewer','',['TEQC Report : ',obs{1},'...'],...
            [600,419,0,-20],'str',str);
    else
        gut('errdlg',['no teqc result : ',rep]);
    end
else
    gut('errdlg',err);
end
set(h,'pointer','arrow');
for n=1:length(obs)
    if ~strcmp(obs{n},obsu{n}), delete(obsu{n}); end
end

% get observation file ---------------------------------------------------------
function [obs,nav]=GetObsFile
data=get(gcf,'userdata'); prm=data.prm;
ts=caltomjd(prm.tstart); te=caltomjd(prm.tend);
rcv=prm.rcv; if length(rcv)>4, rcv=rcv(end-3:end); end
obs={}; nav={};
while ts<=te
    [ts,file,ep]=rinexfile(ts,prm.src,rcv);
    ext={'o','d','o.gz','o.Z','d.gz','d.Z'};
    for n=1:length(ext)
        path=gfilepath(prm.dirs.obs,[file,ext{n}],ep);
        if exist(path), obs={obs{:},path}; break; end
    end
end
ts=caltomjd(prm.tstart);
while ts<=te
    [ts,file,ep]=rinexfile(ts,prm.navsrc,rcv);
    nav={nav{:},gfilepath(prm.dirs.nav,[file,'n'],ep)};
end

function [t,file,ep]=rinexfile(t,src,rcv)
ep=mjdtocal(t); doy=floor(t)-caltomjd([ep(1),1,1])+1;
switch src
case {'rinex','brdc'}
    if strcmp(src,'brdc'), rcv='brdc'; end
    file=sprintf('%s%03d%1d.%02d',rcv,doy,0,mod(ep(1),100));
    t=t+1;
case 'rinex3'
    file=sprintf('%s%03d%1d.%02d',rcv,doy,floor(ep(4)/3)+1,mod(ep(1),100));
    t=t+3/24;
case 'rinex1'
    file=sprintf('%s%03d%c.%02d',rcv,doy,'a'+ep(5),mod(ep(1),100));
    t=t+1/24;
case 'rinexh'
    file=sprintf('%s%03d%c%02d.%02d',rcv,doy,'a'+ep(4),floor(ep(5)/15)*15,mod(ep(1),100));
    t=t+15/1440;
end

% mean without nan -------------------------------------------------------------
function mx=meann(x)
i=find(~isnan(x)); if isempty(i), mx=nan; else mx=mean(x(i)); end

% rms without nan -------------------------------------------------------------
function mx=rmsn(x)
i=find(~isnan(x)); if isempty(i), mx=nan; else mx=sqrt(mean(x(i).^2)); end

% mjd->gpsweek/day/hour --------------------------------------------------------
function wd=gpswd(t)
d=caltomjd(t)-44244; wd(1)=floor(d/7); wd(2)=floor(d-wd(1)*7);

% time string ------------------------------------------------------------------
function str=tstr(td,time)
str=sprintf('%04d/%02d/%02d %02d:%02d:%02.0f',mjdtocal(td,time));
