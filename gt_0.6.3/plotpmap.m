function h=plotpmap(varargin)
%-------------------------------------------------------------------------------
% [system] : GpsTools
% [module] : plot pwv parameters
% [func]   : plot pwv parameters
% [argin]  : none
% [argout] : h = figure handle
% [note]   :
% [version]: $Revision: 4 $ $Date: 06/07/18 14:05 $
%            Copyright(c) 2004-2006 by T.Takasu, all rights reserved
% [history]: 05/11/24  0.1  new
%-------------------------------------------------------------------------------
if nargin>0&strncmp(varargin{1},'cb',2), feval(varargin{:})
else h=PlotMain(varargin{:}); end

% main -------------------------------------------------------------------------
function h=PlotMain(varargin)
h=gut('newfigg','','GPS PWV',[600,400,0,-68],[mfilename,' cbClose']);
if isempty(h), return, end
prm=loadprm('prm_plotpmap','prm_plotpmap_def');
gut('newmenu','data','&Data ',...
    {'&Read...','&Export Plot...','-','&Map Area...','&Options...'},{},...
    {[mfilename,' cbRead'],[mfilename,' cbExport'],'',[mfilename,' cbMap'],...
     [mfilename,' cbOpt']});
gut('newmenu','plot','&Plot',{'&PWV Map'},{'pmap'},[mfilename,' cbPlot']);
gut('newmenu','time','&Time',{},{},[]);
data.time=[]; data.n=[]; data.rcvs={}; data.data=[]; data.gpos=[];
data.prm=prm; set(gcf,'userdata',data);

% callback on close ------------------------------------------------------------
function cbClose
data=get(gcf,'userdata'); prm=data.prm;
[path,file]=fileparts(which(mfilename));
save(fullfile(path,'settings','prm_plotpmap.mat'),'prm');
closereq

% callback on menu read --------------------------------------------------------
function cbRead
data=get(gcf,'userdata'); prm=data.prm;
prm1={
't','Time Data',[]
's','Time Start:End (h)',   [0,0]
'n','Time Interval (sec)',  3600
'b','Stations',[mfilename,' cbSelRcv']
};
prm2={
' ','PWV Data Directory',''
'd','',prm.dirs.est
};
gut('newdlg','','Read Data',[314,185]);
gut('newprms','prm1',[12,73,295,27,118],prm1);
gut('newprms','prm2',[12,33,295,18],prm2);
gut('newokcancelbtn','',[127,4,180,23]);
gut('setprms','prm1',prm.ep,prm.ts,prm.tint'');
gut('setudata','prm1_2',prm.rcvs);
if ~gut('waitok'), return, end
[prm.ep,prm.ts,prm.tint,q]=gut('getprms','prm1');
[q,prm.dirs.est]=gut('getprms','prm2');
prm.rcvs=gut('getudata','prm1_2');
close
data.prm=prm; set(gcf,'userdata',data); ReadData; UpdatePlot;

% callback on menu select station ---------------------------------------------
function cbSelRcv
[rcvs,ok]=editlist('Stations',get(gcbo,'userdata'),'rcv');
if ok, set(gcbo,'userdata',rcvs); end

% callback on menu export -----------------------------------------------------
function cbExport
persistent dirs
if isempty(dirs), dirs=''; end
data=get(gcf,'userdata'); prm=data.prm;
prm1={'','Output Directory', '';'d','',dirs};
gut('newdlg','','Export Plots',[320,76]);
gut('newprms','prm1',[12,35,300,20,120],prm1);
gut('newokcancelbtn','',[133,4,180,23]);
if ~gut('waitok'), return, end
[q,dirs]=gut('getprms','prm1');
close
format={'-djpeg','-r0'}; opts={'-noui'};
f=gcf;
for n=1:length(data.time)
    data.n=n; ep=mjdtocal(caltomjd(prm.ep),data.time(n));
    set(f,'userdata',data); UpdatePlot;
    file=fullfile(dirs,sprintf('pwvmap_%04d%02d%02d%02d%02d.jpg',ep(1:5)));
    print(format{:},opts{:},file);
end

% callback on menu options -----------------------------------------------------
function cbOpt
data=get(gcf,'userdata'); prm=data.prm;
sel1={{'OFF','ON'},0:1};
sel2={{'Right','Bottom'},0:1};
sel3={{'Linear','Cubic','V4'},0:2};
prm1={
'p','Show Station Position',sel1
'p','Show PWV Value',sel1
's','Grid Longitude Min:Max (deg)',[0,0]
's','Grid Latitude Min:Max (deg)',[0,0]
'n','Grid Interval (deg)',0.1
'p','Grid Interpolation',sel3
's','Range Min:Step:Max (mm)',[0,0.5,80]
'p','Color Bar Position',sel2
'c','Station Position Color', 'r'
'c','Contour Line Color', 'none'
};
gut('newdlg','','Options',[320,274]);
gut('newprms','prm1',[12,35,300,23,120],prm1);
gut('newokcancelbtn','',[133,4,180,23]);
gut('setprms','prm1',prm.showf(1),prm.showf(2),prm.area(1:2),prm.area(3:4),...
    prm.gint,prm.intp,prm.range,prm.cbpos,prm.pcolor,prm.ccolor);
if ~gut('waitok'), return, end
[prm.showf(1),prm.showf(2),prm.area(1:2),prm.area(3:4),prm.gint,prm.intp,...
 prm.range,prm.cbpos,prm.pcolor,prm.ccolor]=gut('getprms','prm1');
close
data.prm=prm; set(gcf,'userdata',data); UpdatePlot;

% callback on menu map area ----------------------------------------------------
function cbMap
data=get(gcf,'userdata'); prm=data.prm;
[prm.map,ok]=editmap('Map Area',prm.map);
if ok, data.prm=prm; set(gcf,'userdata',data); UpdatePlot; end

% callback on menu plot --------------------------------------------------------
function cbPlot, UpdatePlot;

% callback on change time ------------------------------------------------------
function cbChgTime
data=get(gcf,'userdata');
data.n=str2num(get(gcbo,'userdata'));
gut('unchkmenu','time')
gut('chkmenu','time',num2str(data.n))
set(gcf,'userdata',data);
UpdatePlot

% plot pwv map -----------------------------------------------------------------
function UpdatePlot
data=get(gcf,'userdata'); prm=data.prm;
pwv=data.data(data.n,:)';
i=find(~isnan(pwv));
pwv=pwv(i,:); gpos=data.gpos(i,:);

clf
[fn,fs]=gut('getfont','g');
p=get(gcf,'position'); tpos=(p(4)-20)/p(4); pos=[0.001,0.001,0.998,tpos];
if prm.cbpos, pos(2)=tpos*0.05; pos(4)=tpos*0.975; else pos(3)=0.92; end

gmt('mmap','proj',prm.map.proj,'cent',prm.map.cent,'base',prm.map.base,...
    'scale',prm.map.scale,'pos',pos,'color',prm.map.color{1},...
    'fontname',fn,'fontsize',fs);
caxis(prm.range([1,end]));
if ~isempty(pwv)
    [lon,lat]=meshgrid(prm.area(3):prm.gint:prm.area(4),prm.area(1):prm.gint:prm.area(2));
    [xl,yl]=gmt('getlim');
    [xs,ys]=gmt('lltoxy',lon,lat);
    [xp,yp]=gmt('lltoxy',gpos(:,2),gpos(:,1));
    switch prm.intp, case 0, intp='cubic'; case 1, intp='linear'; case 2, intp='v4'; end
    zs=griddata(xp,yp,pwv*1E3,xs,ys,intp);
    [c,h]=contourf_v6(xs,ys,zs,prm.range(1):prm.range(2):prm.range(3));
    set(h,'edgecolor',prm.ccolor);
    if ~ischar(prm.map.color{2})|~strcmp(prm.map.color{2},'none');
        gmt('seamask','scolor',prm.map.color{2});
    end
    vpos='middle';
    if prm.showf(1), plot(xp,yp,'.','color',prm.pcolor); vpos='top'; end
    if prm.showf(2)
        for i=find(xl(1)<=xp&xp<=xl(2)&yl(1)<=yp&yp<=yl(2))'
            text(xp(i),yp(i),sprintf('%.1f',pwv(i)*1E3),...
                'horizontal','center','vertical',vpos,'fontname',fn,'fontsize',fs);
        end
    end
end
gmt('mcoast','lcolor','none','scolor','none','ccolor',prm.map.color{3});
gmt('mgrid','gint',prm.map.gint(1),'lint',prm.map.gint(2),'color',prm.map.color{4});
if prm.cbpos, f='colorbarh'; pos=[0.02,tpos*0.04,0.96,0.02];
else f='colorbarv'; pos=[0.95,tpos*0.02,0.015,tpos*0.96]; end
ggt(f,pos,prm.range([1,end]),'','fontname',fn,'fontsize',fs);

title(['GPS-PWV (Precipitable Water Vapor) Map (mm) : ',tstr(prm.ep,data.time(data.n))]);

% read data --------------------------------------------------------------------
function ReadData
d=get(gcf','userdata'); prm=d.prm;
file=fullfile(prm.dirs.est,sprintf('pwvgps_ALL_%04d%02d%02d.mat',prm.ep(1:3)));
if ~exist(file), gt_log(['no pwv estimation file : ',file]); return, end
load(file);
d.time=prm.ts(1)*3600:prm.tint:prm.ts(2)*3600;
[t,i,j]=intersect(time,caltomjd(prm.ep)+d.time/86400);
d.rcvs=rcvs;
d.data=data(i,:);
d.gpos=[];
d.n=1;
poss=readpos(0,0,d.rcvs,'','approx');
for n=1:size(poss,3), d.gpos(n,:)=eceftogeod(poss(1,:,n)'); end

menu={}; val={}; cbs={};
for n=1:length(d.time)
    menu={menu{:},tstr(prm.ep,d.time(n))};
    val={val{:},num2str(n)};
    cbs={cbs{:},[mfilename,' cbChgTime']};
end
delete(gut('geth','time'));
gut('newmenu','time','&Time',menu,val,cbs);
gut('chkmenu','time',num2str(d.n));
set(gcf,'userdata',d,'name',['GPS PWV : ',file])

% time string ------------------------------------------------------------------
function s=tstr(ep,ts)
e=mjdtocal(caltomjd(ep),ts);
s=sprintf('%4d/%02d/%02d %02d:%02d',e(1:5));

% matlab 6 compatible contourf -------------------------------------------------
function [c,h]=contourf_v6(varargin)
v=version; v=str2num(v(1));
if v<=6, [c,h]=contourf(varargin{:}); else [c,h]=contourf('v6',varargin{:}); end
