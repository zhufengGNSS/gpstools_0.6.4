function h=plotion(varargin)
%-------------------------------------------------------------------------------
% [system] : GpsTools
% [module] : plot global tec map
% [func]   : plot global tec map
% [argin]  : 'prm',prm   = parameters
% [argout] : none
% [note]   :
% [version]: $Revision: 16 $ $Date: 2008-12-12 15:49:30 +0900 (金, 12 12 2008) $
%            Copyright(c) 2004-2006 by T.Takasu, all rights reserved
% [history]: 04/12/06  0.1  new
%            05/05/19  0.2  add 3d view
%            08/09/30  0.3  modify apperence (gt_0.6.4)
%-------------------------------------------------------------------------------
if nargin>0&strncmp(varargin{1},'cb',2), feval(varargin{1});
else h=PlotMain(varargin{:}); end

% main -------------------------------------------------------------------------
function h=PlotMain(varargin)
h=gut('newfigg','','Ionospheric Parameters',[600,400,0,-68],[mfilename,' cbClose']);
if isempty(h), return, end
prm=loadprm('prm_plotion','prm_plotion_def');
if nargin>=2&strcmp(varargin{1},'prm'), prm=varargin{2}; end
set(h,'renderermode','auto');
gut('newmenu','data','&Data ',...
    {'&Read...','&New Window','-','&Export Plot...','-','&Map Area','&Options...'},{},...
    {[mfilename,' cbRead'],[mfilename,' cbNew'],'',[mfilename,' cbExport'],'',...
     [mfilename,' cbMap'],[mfilename,' cbOpt']});
gut('newmenu','plot','&Plot',{'&VTEC Map'},{''},{''});
gut('newmenu','time','&Time',{},{},{});
gut('newmenu','height','&Height',{},{},{});
set(gut('newbtnh','sbtn',[0,0,36,16],{'<','>'},[mfilename,' cbSwBtn']),...
    'fontsize',8);
set(gcf,'resizefcn',[mfilename,' cbResize']);
LocatePos;
data.tec=[]; data.prm=prm; set(h,'userdata',data);
if length(varargin)>0, ReadData; UpdatePlot; end

% callback on new --------------------------------------------------------------
function cbNew
data=get(gcf,'userdata');
feval(mfilename,'prm',data.prm);

% callback on close ------------------------------------------------------------
function cbClose
data=get(gcf,'userdata'); prm=data.prm;
[path,file]=fileparts(which(mfilename));
save(fullfile(path,'settings','prm_plotion.mat'),'prm');
closereq

% callback on resize -----------------------------------------------------------
function cbResize, LocatePos

% callback on menu read --------------------------------------------------------
function cbRead
sel1={{'Estimated(Forward)','Estimated(Backward)','Estimated(Smoothed)'},...
      {'ionf','ionb','ionfb'}};
sel2={{'IGS Final','IGS Rapid'},{'igs','igr'}};
data=get(gcf,'userdata'); prm=data.prm;
prm1={
't','Start Time (GPST)',[]
't','End Time (GPST)',  []
'n','Time Interval (sec)', 3600
'n','Processing Unit Time (hr)',  24
'p','Ionos Parameters',{{sel1{1}{:},sel2{1}{:}},{sel1{2}{:},sel2{2}{:}}}
'p','Reference Ionos Parameters',{{''},{''}}
};
prm2={
' ','Ionos Parameters Directory',''
'd','',prm.dirs.est
' ','Reference Ionos Parameters Directory',''
'd','',prm.dirs.ion
};
gut('newdlg','','Read Data',[314,280]);
gut('newprms','prm1',[12,114,295,27,118],prm1);
gut('newprms','prm2',[12,34,295,18],prm2);
gut('newokcancelbtn','',[127,4,180,23]);
gut('setprms','prm1',prm.tstart,prm.tend,prm.tint,prm.tunit,prm.fb,prm.ref);
if ~gut('waitok'), return, end
[prm.tstart,prm.tend,prm.tint,prm.tunit,prm.fb,prm.ref]=gut('getprms','prm1');
[q,prm.dirs.est,q,prm.dirs.ion]=gut('getprms','prm2');
close
data.prm=prm; set(gcf,'userdata',data);
ReadData; UpdatePlot;

% callback on menu export ------------------------------------------------------
function cbExport
persistent file dirs
if isempty(file), file='plotion'; end, if isempty(dirs), dirs=''; end
f=gcf; data=get(f,'userdata');
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
m=data.n;
for n=1:size(data.tec,4)
    data.n=n; set(f,'userdata',data)
    figure(f), UpdatePlot;
    print(format{:},opts{:},sprintf('%s%02d',fs,n));
end
data.n=m; set(f,'userdata',data); UpdatePlot;

% callback on menu map area -----------------------------------------------------
function cbMap
data=get(gcf,'userdata'); prm=data.prm;
[prm.map,ok]=editmap('Map Area',prm.map,[1,1,1,1,1,1,1]);
if ok, data.prm=prm; set(gcf,'userdata',data); UpdatePlot; end

% callback on menu options ------------------------------------------------------
function cbOpt
data=get(gcf,'userdata'); prm=data.prm;
sel1={{'OFF','ON'},0:1};
sel2={{'Contour Line','Contour Filled','Surface'},0:2};
prm1={
'p','Plot Type',               sel2
'c','Contour Color',           ''
'c','Geomagnetic Equator Color',''
'p','Show Contour Label',      sel1
's','Range Min:Step:Max',      [0,2.5,100]
'p','Fix Local Time',          sel1
's','Light Angle Az:El (deg)', [0,0]
's','Diffuse/Ambient Light Str.',[0,0]
};
q='';
gut('newdlg','','Options',[280,225]);
gut('newprms','prm1',[15,35,260,23,120],prm1);
gut('newokcancelbtn','',[100,4,170,23]);
gut('setprms','prm1',prm.ptype,prm.ccolor,prm.ecolor,prm.popt(2),prm.range,...
    prm.fixlt,prm.light,prm.lstr);
if ~gut('waitok'), return, end
[prm.ptype,prm.ccolor,prm.ecolor,prm.popt(2),prm.range,prm.fixlt,prm.light,...
 prm.lstr]=gut('getprms','prm1');
close
data.prm=prm; set(gcf,'userdata',data);
UpdatePlot;

% callback on change time ------------------------------------------------------
function cbChgTime
data=get(gcf,'userdata');
data.n=str2num(get(gcbo,'userdata'));
gut('unchkmenu','time')
gut('chkmenu','time',num2str(data.n))
set(gcf,'userdata',data);
UpdatePlot

% callback on change height ----------------------------------------------------
function cbChgHeight
data=get(gcf,'userdata');
data.m=str2num(get(gcbo,'userdata'));
gut('unchkmenu','height')
gut('chkmenu','height',num2str(data.m))
set(gcf,'userdata',data);
UpdatePlot

% callback on push swithch button ----------------------------------------------
function cbSwBtn
data=get(gcf,'userdata');
n=size(data.tec,4);
switch get(gcbo,'tag')
case 'sbtn_1', if data.n>1, data.n=data.n-1; end
case 'sbtn_2', if data.n<n, data.n=data.n+1; end
end
gut('unchkmenu','time')
gut('chkmenu','time',num2str(data.n))
set(gcf,'userdata',data);
UpdatePlot;

% update plot ------------------------------------------------------------------
function UpdatePlot
data=get(gcf,'userdata'); prm=data.prm;
clf, axis off
if isempty(data.tec), return, end
t=data.time(data.n)/3600; cent=prm.map.cent; base=prm.map.base;
if prm.fixlt, cent(1)=cent(1)-t*15; base(1)=base(1)-t*15; end

if prm.ptype==2, m=0.08; else m=0; end
p=get(gcf,'position'); tpos=(p(4)-20)/p(4);pos=[0.006,0.008,0.988-m,tpos-0.008];
%p=get(gcf,'position'); tpos=(p(4)-20)/p(4);pos=[0,0,1,1];
[fn,fs]=gut('getfont','g');

gmt('mmap','proj',prm.map.proj,'cent',cent,'base',base,'scale',prm.map.scale,...
    'pos',pos,'fontname',fn,'fontsize',fs);
h=gmt('mcoast','lcolor',prm.map.color{1},'scolor',prm.map.color{2},'ccolor',...
       prm.map.color{3});
h=[h;gmt('mgrid','gint',prm.map.gint(1),'lint',prm.map.gint(2),'color',prm.map.color{4})];

vs=data.tec(:,2:end,data.m,data.n);
[lon,lat]=meshgrid(data.lons(2:end),data.lats);
[xs,ys,zs]=gmt('lltoxy',lon,lat);
switch prm.map.proj
case {'eq','miller','mercat','mollw'}
    [x,i]=sort(xs(1,:),2); xs=xs(:,i); ys=ys(:,i); vs=vs(:,i);
    dx=xs(:,2)-xs(:,1); xs=[xs(:,1)-dx,xs];
    ys=[ys(:,end),ys]; vs=[vs(:,end),vs]; zs=vs;
case {'ortho','polars','polarn','lambert'}
    [xs,ys,zs]=gmt('lltoxy',lon,lat);
    [x,i]=sort(gmt('normlon',lon(1,:)-cent(1)),2); i=[i(end),i];
    xs=xs(:,i); ys=ys(:,i); vs=vs(:,i); zs=zs(:,i);
    if strcmp(prm.map.proj,'ortho'), vs(zs<0)=nan; end
end
if ~strcmp(prm.ccolor,'none'), lc=prm.ccolor; else lc='k'; end
switch prm.ptype
case 0
    [c,h]=contour_v6(xs,ys,vs,prm.range(1):prm.range(2):prm.range(3));
    if ~strcmp(prm.ccolor,'none')
        set(h,'edgecolor',prm.ccolor);
    end
    if prm.popt(2)
        clabel(c,h,'fontname',fn,'fontsize',fs,'color',lc);
    end
case 1
    [c,h]=contourf_v6(xs,ys,vs,prm.range(1):prm.range(2):prm.range(3));
    set(h,'edgecolor',prm.ccolor);
    if prm.popt(2)
        clabel(c,h,'fontname',fn,'fontsize',fs,'color',lc);
    end
case 2
    surf(xs,ys,zs,vs,'edgecolor','none','facecolor','interp','facelighting','phong',...
         'diffusestrength',prm.lstr(1),'ambientstrength',prm.lstr(2),...
         'specularstrength',0.5,'specularexponent',1,'specularcolorreflectance',0.5);
    lightangle(prm.light(1),prm.light(2));
end

if ~strcmp(prm.ecolor,'none')
    lats=[]; lons=[];
    for lonm=0:5:360
        [lat,lon]=geomtogeoc(0,lonm*pi/180); lats=[lats;lat]; lons=[lons;lon];
    end
    h=[h;gmt('mplot',lons*180/pi,lats*180/pi,prm.ecolor)];
end
if strcmp(prm.map.proj,'ortho'), gmt('setz',h,1); else gmt('setz',h,prm.range(3)); end
ti=sprintf('Ionosphere Map VTEC(TECU) : %s H=%.0fkm',...
           tstr(mjdtocal(data.td,data.time(data.n))),data.hgts(data.m));
title(ti);

caxis(prm.range([1,3]))
if prm.ptype==2
    ggt('colorbarv',[0.94,0.015,0.015,tpos-0.03],prm.range([1,3]),'',...
        'fontname',fn,'fontsize',fs);
end

% read data --------------------------------------------------------------------
function ReadData
data=ReadIonFiles(get(gcf,'userdata')); prm=data.prm;
menu1={}; val1={}; cbs1={};
for n=1:length(data.time)
    menu1={menu1{:},tstr(mjdtocal(data.td,data.time(n)))};
    val1={val1{:},num2str(n)};
    cbs1={cbs1{:},[mfilename,' cbChgTime']};
end
menu2={}; val2={}; cbs2={};
for n=1:length(data.hgts)
    menu2={menu2{:},['H=',num2str(data.hgts(n)),'km']};
    val2={val2{:},num2str(n)};
    cbs2={cbs2{:},[mfilename,' cbChgHeight']};
end
delete(gut('geth','time'));
delete(gut('geth','height'));
gut('newmenu','time','&Time',menu1,val1,cbs1);
gut('newmenu','height','&Height',menu2,val2,cbs2);
gut('chkmenu','time',num2str(data.n));
gut('chkmenu','height',num2str(data.m));
ti=['Ionospheric Parameters : ',prm.dirs.est,' (',prm.fb,')'];
set(gcf,'userdata',data,'name',ti);

% read ionosphere files --------------------------------------------------------
function data=ReadIonFiles(data)
prm=data.prm;
[data.td,ts]=caltomjd(prm.tstart); [tn,te]=caltomjd(prm.tend);
te=te+(tn-data.td)*86400;
h=gcf; set(h,'pointer','watch')
data.time=[]; data.tec=[]; data.lats=[]; data.lons=[]; data.hgts=[];
data.n=[]; data.m=[];
for td=data.td:tn
    ep=mjdtocal(td); doy=td-caltomjd([ep(1),1,1])+1;
    switch prm.fb
    case 'igs', file=sprintf('igsg%03d0.%02di',doy,mod(ep(1),100));
    case 'igr', file=sprintf('igrg%03d0.%02di',doy,mod(ep(1),100));
    otherwise,  file=sprintf('%s_%04d%02d%02d%02d.mat',prm.fb,ep(1:3),0);
    end
    file=gfilepath(prm.dirs.est,file,ep);
    [epoch,time,tec,rms,lats,lons,hgts]=readionex(file);
    siz=size(tec);
    if isempty(epoch)
        ;
    elseif siz(1)>0&isempty(data.tec)
        [tdd,ts]=caltomjd(epoch);
        data.time=ts+time+(tdd-data.td)*86400;
        data.tec=tec;
        data.lats=lats;
        data.lons=lons;
        data.hgts=hgts;
        sizi=size(tec); 
    elseif all(siz(1:3)==sizi(1:3))
        [tdd,ts]=caltomjd(epoch);
        data.time=[data.time;time+ts+(tdd-data.td)*86400];
        for n=1:size(data.tec,3)
            data.tec(:,:,n,1:end+siz(4))=cat(4,data.tec(:,:,n,:),tec(:,:,n,:));
        end
    end
end
set(h,'pointer','arrow')
if isempty(data.tec), return, end
i=find(ts<=data.time&data.time<=te&mod(data.time,prm.tint)==0);
[t,j]=unique(data.time(i));
data.time=data.time(i(j));
data.tec=data.tec(:,:,:,i(j));
data.n=1;
data.m=1;

% locate gui objects -----------------------------------------------------------
function LocatePos
p=get(gcf,'position'); m=15;
gut('setpos','sbtn_1',[p(3)-m-32,p(4)-16]);
gut('setpos','sbtn_2',[p(3)-m-16,p(4)-16]);

% time string ------------------------------------------------------------------
function str=tstr(ep), str=sprintf('%04d/%02d/%02d %02d:%02d:%02.0f',ep);

% matlab 6 compatible plots ----------------------------------------------------
function [c,h]=contour_v6(varargin)
v=version; v=str2num(v(1));
if v<=6, [c,h]=contour(varargin{:}); else [c,h]=contour('v6',varargin{:}); end

function [c,h]=contourf_v6(varargin)
v=version; v=str2num(v(1));
if v<=6, [c,h]=contourf(varargin{:}); else [c,h]=contourf('v6',varargin{:}); end
