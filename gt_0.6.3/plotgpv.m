function h=plotgpv(varargin)
%-------------------------------------------------------------------------------
% [system] : GpsTools
% [module] : plot gpv data
% [func]   : plot gpv data
% [argin]  : 'prm',prm   : parameters
% [argout] : h = figure handle
% [note]   :
% [version]: $Revision: 16 $ $Date: 06/07/19 5:26 $
%            Copyright(c) 2004-2006 by T.Takasu, all rights reserved
% [history]: 05/05/06  0.1  new
%-------------------------------------------------------------------------------
if nargin>0&strncmp(varargin{1},'cb',2), feval(varargin{1});
else h=PlotMain(varargin{:}); end

% main -------------------------------------------------------------------------
function h=PlotMain(varargin)
h=gut('newfigg','','Grid Point Values',[600,400,0,-68],[mfilename,' cbClose']);
if isempty(h), return, end
set(h,'render','painter','renderermode','auto');
prm=loadprm('prm_plotgpv','prm_plotgpv_def');
if nargin>=2&strcmp(varargin{1},'prm'), prm=varargin{2}; end
gut('newmenu','data','&Data ',...
    {'&Read...','&New Window','-','&Export Plot...','-','&Map Area...',...
     '&Color Map...','&Options...'},{},...
    {[mfilename,' cbRead'],[mfilename,' cbNew'],'',[mfilename,' cbExport'],'',...
     [mfilename,' cbMap'],[mfilename,' cbCmap'],[mfilename,' cbOpt']});
gut('newmenu','plot','&Plot',...
    {'&Map','-','&Contour Line','Contour &Filled','&Mesh Grid','&Surface',...
     '&Vector Arrow','Stream &Line','-','Contour Line 2','Contour Filled 2',...
      'Vector Arrow 2','Stream Line 2','-','Color &Bar'},...
    {'map','','contl1','contf1','mesh1','surf1','vect1','stream1','','contl2',...
     'contf2','vect2','stream2','','cbar'},[mfilename,' cbPlot']);
gut('newmenu','height','&Height',{},{},'');
gut('newmenu','time','&Time',{},{},'');
gut('chkmenu','plot',prm.ptype);
set(gut('newtext','pos',[0,0,180,15],'',2),'visible','off','backgroundcolor','w');
set(gut('newbtnh','tbtn',[0,0,69,16],{'.','<','>','.'},[mfilename,' cbToolBtn']),...
    'fontsize',8);
set(gcf,'resizefcn',[mfilename,' cbResize'],'windowbuttonmotionfcn',[mfilename,' cbBtnMove']);
LocatePos;
data.time=[]; data.ft=[]; data.data1={}; data.data2={}; data.height1=[];
data.height2=[]; data.gprm=[]; data.mplot=[]; data.gplot={}; data.cplot=[];
data.ti={'',''}; data.cs={}; data.lock=0;
data.prm=prm; set(gcf,'userdata',data);
if ~isempty(prm.cmap), set(h,'colormap',prm.cmap); end
if length(varargin)>0, ReadData; UpdatePlot; end

% callback on new --------------------------------------------------------------
function cbNew
data=get(gcf,'userdata');
feval(mfilename,'prm',data.prm);

% callback on close ------------------------------------------------------------
function cbClose
data=get(gcf,'userdata'); prm=data.prm;
[path,file]=fileparts(which(mfilename));
save(fullfile(path,'settings','prm_plotgpv.mat'),'prm');
closereq

% callback on resize -----------------------------------------------------------
function cbResize, LocatePos

% callback on menu read --------------------------------------------------------
function cbRead
sel1={{'','Pressure MSL','Geopotential Height','Temperature','Wind',...
       'Vertical Velocity','Relative Humidity','Rain','Cloud','Geoid Height'},...
      {'','pmsl','geoh','temp','wind','vvel','humi','rain','clou','geoid'}};
sel2={{'JMA MSM','JMA RSM','JMA GSM','JMA MSM Online','JMA RSM Online',...
       'JMA GSM Online','GSI','EGM96'},...
      {'msm','rsm','gsm','mso','rso','gso','gsi','egm'}};
sel3={{'All','Surface','10hPa','20hPa','30hPa','50hPa','70hPa','100hPa',...
       '150hPa','200hPa','250hPa','300hPa','400hPa','500hPa','600hPa','700hPa',...
       '800hPa','850hPa','900hPa','925hPa','950hPa','975hPa','1000hPa'},...
      [nan,0,10,20,30,50,70,100,150,200,250,300,400,500,600,700,800,850,900,...
       925,950,975,1000]};
data=get(gcf,'userdata'); prm=data.prm;
prm1={
't','Start Time (UTC)',  []
't','End   Time (UTC)',  []
'n','Initial Time Interval (hr)', 6
'n','Forecast Time Interval (hr)', 0
'p','Gpv Data Model',     sel2
'p','Gpv Data Type 1',    sel1
'p','Gpv Data Type 2',    sel1
'p','Gpv Data Height 1',  sel3
'p','Gpv Data Height 2',  sel3
};
prm2={
' ','Gpv Data Directory',''
'd','',prm.dirs.gpv
};
gut('newdlg','','Read Data',[314,308]);
gut('newprms','prm1',[12,78,295,25,118],prm1);
gut('newprms','prm2',[12,34,295,18,118],prm2);
gut('newokcancelbtn','',[127,4,180,23]);
gut('setprms','prm1',prm.tstart,prm.tend,prm.tint,prm.ft,prm.src,prm.type{1},...
    prm.type{2},prm.height(1),prm.height(2));
if ~gut('waitok'), return, end
[prm.tstart,prm.tend,prm.tint,prm.ft,prm.src,prm.type{1},prm.type{2},...
 prm.height(1),prm.height(2)]=gut('getprms','prm1');
[q,prm.dirs.gpv]=gut('getprms','prm2');
close
data.prm=prm; set(gcf,'userdata',data); ReadData; UpdatePlot;

% callback on menu export ------------------------------------------------------
function cbExport
persistent file dirs
if isempty(file), file='plotgpv'; end, if isempty(dirs), dirs=''; end
sel1={{'Single Frame','All Time Frames'},0:1};
sel2={{'EPS (*.eps)','JPEG (*.jpg)','TIFF (*.tif)','META FILE (*.emf)'},...
      {{'-depsc','-tiff'},{'-djpeg','-r0'},{'-dtiff','-r0'},{'-dmeta'}}};
prm1={
'p','Frame Sequence',sel1
'e','File Name',file
'p','File Format',sel2
' ','Directory','';'d','',dirs
};
gut('newdlg','','Export Plots',[274,160]);
gut('newprms','prm1',[12,34,255,24,150],prm1);
gut('newokcancelbtn','',[88,4,180,23]);
if ~gut('waitok'), return, end
[frame,file,format,q,dirs]=gut('getprms','prm1');
close
fs=fullfile(dirs,file); opts={'-noui'};
if ~frame
    print(format{:},opts{:},fs);
    return
end
data=get(gcf,'userdata');
for n=1:length(data.time)
    gut('unchkmenu','time');
    gut('chkmenu','time',num2str(n));
    UpdatePlot([0,0,1]);
    print(format{:},opts{:},sprintf('%s%02d',fs,n));
end

% callback on menu options -----------------------------------------------------
function cbOpt
sel1={{'OFF','ON'},0:1};
sel2={{'Right','Bottom'},0:1};
sel3={{'Gpv Data 1','Gpv Data 2'},1:2};
sel4={{'Forward','Backward','Both'},0:2};
sel5={{'Z Axis Min','Z Axis Max','On Surface/Mesh'},1:3};
sel6={{'Flat','Adjust To Mesh/Surface'},0:1};
[az,el]=view; data=get(gcf,'userdata'); prm=data.prm;
prm1={
'p','Perspective Projection',  sel1;
's','View Angle Az:El (deg)',  [0,0]
'p','Show Axis',               sel1
'p','Show Box',                sel1
's','Z Factor:Min:Max (inf:auto)',[0,0,0]
'p','Z Position of Map',       sel5
'p','Z Position of Data',      sel5
'p','Color Bar Postion',       sel2
'p','Color Bar Data',          sel3
'c','Contour 1 Line Color',    ''
'p','Contour 1 Line Label',    sel1
's','Cont 1 Min:Step:Max (0:auto)',[0,0,0]
'c','Contour 2 Line Color',         ''
'p','Contour 2 Line Label',         sel1
's','Cont 2 Min:Step:Max (0:auto)',[0,0,0]
};
prm2={
'c','Background Color',        ''
'c','Mesh Grid Color',         ''
'n','Mesh Grid Interval',      2
'c','Surface Color',           ''
's','Surface Refrect. Parametes',[0,0,0]
'c','Vector Arrow Color',      ''
'n','Vector Arrow Length',     10
'n','Vector Arrow Interval',   5
'c','Stream Line Color',       nan
'n','Stream Line Length',      1000
'p','Stream Line Direction',   sel4
'n','Stream Line Interval',    5
's','Light Angle Az:El (deg)', [0,0]
's','Diffuse/Ambient Light Str.',[0,0]
};
q='';
gut('newdlg','','Options',[558,362]);
gut('newprms','prm1',[15,11,260,23,110],prm1);
gut('newprms','prm2',[290,34,260,23,110],prm2);
gut('newokcancelbtn','',[380,4,170,23]);
gut('setprms','prm1',prm.vproj,[az,el],prm.axis,prm.box,prm.zlim,prm.zpos(1),...
    prm.zpos(2),prm.cbpos,prm.cbdata,prm.ccolor{1},prm.clabel{1},...
    prm.cstep{1},prm.ccolor{2},prm.clabel{2},prm.cstep{2});
gut('setprms','prm2',prm.bcolor,prm.mcolor,prm.mint,prm.fcolor,prm.fparam,...
    prm.vcolor,prm.vlen,prm.vint,prm.scolor,prm.slen,prm.sdir,prm.sint,...
    prm.light,prm.lstr);
if ~gut('waitok'), return, end
[prm.vproj,prm.view,prm.axis,prm.box,prm.zlim,prm.zpos(1),prm.zpos(2),prm.cbpos,...
 prm.cbdata,prm.ccolor{1},prm.clabel{1},prm.cstep{1},prm.ccolor{2},...
 prm.clabel{2},prm.cstep{2}]=gut('getprms','prm1');
[prm.bcolor,prm.mcolor,prm.mint,prm.fcolor,prm.fparam,prm.vcolor,prm.vlen,...
 prm.vint,prm.scolor,prm.slen,prm.sdir,prm.sint,prm.light,prm.lstr]=gut('getprms','prm2');
close
data.prm=prm; set(gcf,'userdata',data);
UpdatePlot;

% callback on menu map area -----------------------------------------------------
function cbMap
data=get(gcf,'userdata'); prm=data.prm;
[prm.map,ok]=editmap('Map Area',prm.map,[1,1,1,1,1,1,1]);
if ok, data.prm=prm; set(gcf,'userdata',data); UpdatePlot([1,1,0]); end

% callback on menu color map ----------------------------------------------------
function cbCmap
data=get(gcf,'userdata'); prm=data.prm;
[ok,prm.cmap]=ggt('cmapdlg','','Color Map',get(gcf,'colormap'));
if ~ok, return, end
data.prm=prm; set(gcf,'userdata',data,'colormap',prm.cmap);
UpdatePlot([0,1,0]);

% callback on menu plot ---------------------------------------------------------
function cbPlot
data=get(gcf,'userdata'); prm=data.prm;
gut('togglechk',gcbo);
prm.ptype=gut('getchk','plot');
data.prm=prm; set(gcf,'userdata',data);
UpdatePlot([1,1,0])

% callback on change time -------------------------------------------------------
function cbTime
persistent lock brk
switch get(gcbo,'userdata')
case 'movie'
    if ~isempty(lock), brk=1; return, end
    lock=1; brk=0;
    data=get(gcf,'userdata');
    for n=1:length(data.time)
        gut('unchkmenu','time');
        gut('chkmenu','time',num2str(n));
        UpdatePlot([0,0,1]); drawnow;
        if brk, break, end
    end
    lock=[];
otherwise
    gut('unchkmenu','time');
    set(gcbo,'checked','on'); 
    UpdatePlot([0,0,1])
end
UpdateTimeF

% callback on tool button ------------------------------------------------------
function cbToolBtn
data=get(gcf,'userdata');
n=length(data.time); if n<1, return, end
t=gut('getchk','time'); m=str2num(t{1}); 
switch get(gcbo,'tag')
case 'tbtn_1', m=1;
case 'tbtn_2', m=max(1,m-1);
case 'tbtn_3', m=min(n,m+1);
case 'tbtn_4', m=n;
end
gut('unchkmenu','time');
gut('chkmenu','time',num2str(m));
UpdateTimeF
UpdatePlot([0,0,1]);

% update forward time ----------------------------------------------------------
function UpdateTimeF
data=get(gcf,'userdata'); n=length(data.time);
t=gut('getchk','time'); m=str2num(t{1});
if m<=1, gut('setdis',{'tbtn_1','tbtn_2'}); else gut('setena',{'tbtn_1','tbtn_2'}); end
if m>=n, gut('setdis',{'tbtn_3','tbtn_4'}); else gut('setena',{'tbtn_3','tbtn_4'}); end

% callback on change height ----------------------------------------------------
function cbHeight
height=get(gcbo,'userdata');
for h=get(get(gcbo,'parent'),'children')'
    if strncmp(get(h,'userdata'),height,2), set(h,'checked','off'); end
end
set(gcbo,'checked','on'); 
UpdatePlot([1,1,0])

% update plot ------------------------------------------------------------------
function UpdatePlot(update)
if nargin<1, update=[1,1,1]; end
data=get(gcf,'userdata'); prm=data.prm;
set(gcf,'pointer','watch');
if isempty(data.data1), return, end
p=get(gcf,'position'); tpos=(p(4)-20)/p(4); pos=[0.001,0.001,0.998,tpos];
[fn,fs]=gut('getfont','g'); 
if any(strcmp(prm.ptype,'cbar'))
    if prm.cbpos, pos(2)=tpos*0.05; pos(4)=tpos*0.975; else pos(3)=0.92; end
end
if prm.zlim(1)>0, zfact=prm.zlim(1); else zfact=0.1; end
gmt('mmap','proj',prm.map.proj,'cent',prm.map.cent,'base',prm.map.base,...
    'scale',prm.map.scale,'pos',pos,'zfact',zfact,'fontname',fn,...
    'fontsize',fs,'color',prm.bcolor);
if prm.axis, axis on, else axis off; end, if prm.box, box on; else box off; end
h=gut('getchk','time'); if ~isempty(h), t=str2num(h{1}); else t=[]; end
if any(update(1:2))
    for n=1:length(data.gplot), delete(data.gplot{n}); end
    data.gplot={}; for n=1:length(data.time), data.gplot{n}=[]; end
end
if update(3) % update time
    for n=1:length(data.time), set(data.gplot{n},'visible','off'); end
    if isempty(data.gplot{t}), update(2)=1;
    else set(data.gplot{t},'visible','on'); end
end
zg=[];
if update(2)&~isempty(data.gprm) % update data
    delete(data.gplot{t}); h=[]; ti=''; he=''; data.cs={[0,1],[0,1]};
    data1=data.data1{t}; data2=data.data2{t}; hh=gut('getchk','height');
    if ~isempty(data1)&length(hh)>=1
        he=hh{1}(3:end); n=find(strcmp(data.height1,he));
        [vs1,vec1,ti,data.cs{1}]=GridData(data1(:,:,n,:),prm.type{1},prm.cstep{1});
        [h1,xg,yg,zg]=PlotData(vs1,vec1,data.cs{1},prm,data.gprm,1);
        h=[h;h1];
        if ~strcmp(prm.map.proj,'ortho')
            zl=[min(min(vs1)),max(max(vs1))]; k=find(~isinf(prm.zlim(2:3))); zl(k)=prm.zlim(k+1);
            if zl(1)<zl(2), zlim(zl); end
        else zlim([-1,1]); end
    end
    if ~isempty(data2)&length(hh)>=2
        he2=hh{2}(3:end); n=find(strcmp(data.height2,he2));
        [vs2,vec2,ti2,data.cs{2}]=GridData(data2(:,:,n,:),prm.type{2},prm.cstep{2});
        h=[h;PlotData(vs2,vec2,data.cs{2},prm,data.gprm,2)];
        if ~strcmp(ti,ti2), ti=[ti,':',ti2]; end
        if ~strcmp(he,he2), he=[he,':',he2]; end
    end
    data.gplot{t}=h; data.ti{1}=ti; data.ti{2}=['H=',he];
end
if update(1) % update map
    delete(data.mplot);
    zl=zlim; if prm.zpos(1)<=2, zp=zl(prm.zpos(1)); else zp=zl(1); end
    if strcmp(prm.map.proj,'ortho')&strcmp(prm.map.color{1},'none')&...
       ~strcmp(prm.map.color{3},'none')&prm.zpos(1)==3&~isempty(zg)
        [lon,lat]=gmt('getcoast'); h=[];
        for n=1:length(lon)
            h=[h;plotonsurf(lon{n},lat{n},xg,yg,zg,data.gprm,prm.map.color{3})];
        end
        for lon=-180:prm.map.gint(1):180
            lat=-90:5:90;
            h=[h;plotonsurf(repmat(lon,size(lat),1),lat,xg,yg,zg,data.gprm,...
                            prm.map.color{4},'linestyle',':')];
        end
        for lat=-90:prm.map.gint(1):90
            lon=-180:5:180; if lat==0, s='-'; else s=':'; end
            h=[h;plotonsurf(lon,repmat(lat,size(lon),1),xg,yg,zg,data.gprm,...
                            prm.map.color{4},'linestyle',s)];
        end
    else
        h=[gmt('mcoast','lcolor',prm.map.color{1},'scolor',prm.map.color{2},...
               'ccolor',prm.map.color{3},'zpos',zp,'ambientstrength',1);
           gmt('mgrid','gint',prm.map.gint(1),'lint',prm.map.gint(2),'color',...
               prm.map.color{4},'zpos',zp)];
    end
    data.mplot=h;
end
if any(update(1:2))
    if any(strcmp(prm.ptype,'map')), on='on'; else on='off'; end
    set(data.mplot,'visible',on);
end
ti=data.ti{1};
if isempty(findstr(data.ti{1},'Geoid'))
    ti=[ti,sprintf(' : %s %s',tstr(prm.td,data.time(t),data.ft(t)),data.ti{2})];
end
title(ti);
if strcmp(prm.map.proj,'ortho')
    view(0,90); camproj('orthographic');
else
    view(prm.view(1),prm.view(2));
    if prm.vproj, camproj('perspective'); else camproj('orthographic'); end
end
if ~isempty(data.cs)
    delete(data.cplot); data.cplot=[];
    if any(strcmp(prm.ptype,'cbar'))&~isempty(data.cs{prm.cbdata})
        if prm.cbpos, f='colorbarh'; pos=[0.02,tpos*0.04,0.96,0.02];
        else f='colorbarv'; pos=[0.95,tpos*0.02,0.015,tpos*0.96]; end
        data.cplot=ggt(f,pos,data.cs{prm.cbdata}([1,end]),'','fontname',fn,'fontsize',fs);
    end
end
set(gcf,'userdata',data,'pointer','arrow');
UpdateTimeF;

% plot on surface --------------------------------------------------------------
function h=plotonsurf(lon,lat,xg,yg,zg,gprm,color,varargin)
[xs,ys]=gmt('lltogrid',lon,lat,gprm);
[xs,ys,zs]=gmt('lltoxy',lon,lat,interp2(xg,yg,zg,xs,ys));
zs(zs<0)=nan; h=plot3(xs,ys,zs,'color',color,varargin{:});

% get grid data ----------------------------------------------------------------
function [vs,vec,ti,cs]=GridData(data,type,cs)
vs=double(data); vec=[];
switch type
case 'pmsl', ti='Pressure MSL(hPa)';    cstep=[980,1,1040];
case 'geoh', ti='Geop Height(m)';       cstep=[0,100,18000];
case 'temp', ti='Temperature(C)';       cstep=[-80,1,40];
case 'vvel', ti='Vert Verocity(hPa/s)'; cstep=[-50,5,50];
case 'humi', ti='Rel Humidity(%)';      cstep=[0,5,100];
case 'rain', ti='Rain(mm/s)';           cstep=[0,5,100];
case 'clou', ti='Cloud';                cstep=[0,50,1000];
case 'geoid',ti='Geoid Height(m)';      cstep=[-100,5,100];
case 'wind'
    ti='Wind (m/sec)'; cstep=[0,2,100];
    vec=squeeze(vs); vs=sqrt(vec(:,:,1).^2+vec(:,:,2).^2);
end
i=find(cs==0); cs(i)=cstep(i);
if cs(1)>=cs(3)|cs(2)<=0, cs=cstep; disp('warning : contour step error'),  end
cs=(floor(cs(1)/cs(2)):ceil(cs(3)/cs(2)))*cs(2);

% plot grid data ---------------------------------------------------------------
function [h,xg,yg,zg]=PlotData(vs,vec,cs,prm,gprm,n)
i=1:gprm.nx; j=1:gprm.ny;
[xg,yg]=meshgrid(i,j);
zg=(vs-cs(end))/(cs(end)-cs(1))*prm.zlim(1); zgi=zg;
[lon,lat]=gmt('gridtoll',xg,yg,gprm);
switch prm.map.proj
case {'ortho','mollw'}
    [xs,ys,zs]=gmt('lltoxy',lon,lat,zg);
    if gprm.type==0
        [x,i]=sort(gmt('normlon',lon(1,:)-prm.map.cent(1)),2);
        if gmt('normlon',gprm.lon1-gprm.dx)==gprm.lon2, i=[i,i(1)]; end
        xs=xs(:,i); ys=ys(:,i); zs=zs(:,i); vs=vs(:,i);
        if ~isempty(vec), vec=vec(:,i,:); end
        lon=lon(:,i); lat=lat(:,i); zgi=zg(:,i);
        gprm.lon1=gprm.lon1+(i(1)-1)*gprm.dx;
    end
    if strcmp(prm.map.proj,'mollw'), zs=vs; end
case 'lambert'
    [xs,ys,zm]=gmt('lltoxy',lon,lat);
    if gprm.type==0
        i=find(zm(2,:)>=0); [x,k]=sort(xs(2,i),2); i=i(k); j=find(lat(:,1)>-60);
        xs=xs(j,i); ys=ys(j,i); vs=vs(j,i);
        if ~isempty(vec), vec=vec(j,i,:); end, lon=lon(j,i); lat=lat(j,i);
        gprm.lon1=gprm.lon1+(i(1)-1)*gprm.dx;
    end
    zs=vs;
otherwise
    [xs,ys]=gmt('lltoxy',lon,lat);
    if gprm.type==0
        [x,i]=sort(xs(1,:),2);
        xs=xs(:,i); ys=ys(:,i); vs=vs(:,i);
        if ~isempty(vec), vec=vec(:,i,:); end, lon=lon(:,i); lat=lat(:,i);
        gprm.lon1=gprm.lon1+(i(1)-1)*gprm.dx;
    end
    zs=vs;
end
h=[]; pt={}; cont=[0,0];
for m=1:length(prm.ptype)
    if ~isempty(prm.ptype{m})&str2num(prm.ptype{m}(end))==n
        pt={pt{:},prm.ptype{m}(1:end-1)};
    end
end
if any(strcmp(pt,'mesh')),  h=[h;PlotMeshGrid(xs,ys,zs,vs,cs,prm)]; end
if any(strcmp(pt,'surf')),  h=[h;PlotSurface(xs,ys,zs,vs,cs,prm)]; end
if any(strcmp(pt,'contf'))&~any(strcmp(pt,'surf')), cont(1)=1; end
if any(strcmp(pt,'contl')), cont(2)=1; end
if any(cont)
    switch prm.map.proj
    case {'ortho','mollw'}
        h=[h;PlotContMap(xs,ys,zs,vs,cs,prm,gprm,cont,n)];
    otherwise
        h=[h;PlotContour(xs,ys,vs,cs,prm,cont,n)];
    end
end
if any(strcmp(pt,'stream'))&~isempty(vec)
     h=[h;PlotStreamLine(vec,prm,gprm)];
end
if any(strcmp(pt,'vect'))&~isempty(vec)
    [xg,yg]=meshgrid(i,j);
    [u,v]=gmt('gridtouv',xg,yg,vec(:,:,1),vec(:,:,2),gprm);
    [xe,ye]=gmt('lltoxy',lon+0.01,lat,zgi); xe=xe-xs; ye=ye-ys; % east
    [xn,yn]=gmt('lltoxy',lon,lat+0.01,zgi); xn=xn-xs; yn=yn-ys; % north
    h=[h;PlotVecArrow(xs,xe,xn,ys,ye,yn,zs,vs,u,v,cs,prm)];
end

% plot contour line/filled -----------------------------------------------------
function h=PlotContour(xs,ys,vs,cs,prm,cont,n)
if cont(1), [c,h]=contourf_v6(xs,ys,vs,cs); else [c,h]=contour_v6(xs,ys,vs,cs); end
if cont(2), cc=prm.ccolor{n}; else cc='none'; end
cm=get(gcf,'colormap'); nc=size(cm,1); s=(cs(end)-cs(1))/(nc-1); zl=zlim;
for hh=h'
    v=get(hh,'userdata');
    m=floor((v-cs(1))/s)+1; if m<1, m=1; elseif m>nc, m=nc; end
    if cont(1), fc=cm(m,:); else fc='none'; end
    if cont(2)&strcmp(cc,'none'), ccc=cm(m,:); else ccc=cc; end
    set(hh,'facecolor',fc,'edgecolor',ccc);
    if prm.zpos(2)<=2, gmt('setz',hh,zl(prm.zpos(2))); else gmt('setz',hh,v); end
end
if cont(2)&prm.clabel{n}&~isempty(c)
    [fn,fs]=gut('getfont','g'); if strcmp(cc,'none'), cc='k'; end
    h=[h;clabel(c,h,'fontname',fn,'fontsize',fs,'color',cc)];
end

% plot contour line/filled in map ----------------------------------------------
function h=PlotContMap(xs,ys,zs,vs,cs,prm,gprm,cont,n)
[c,hs]=contour_v6(vs,cs); hl=[];
if cont(2)&prm.clabel{n}&~isempty(c)
    hl=clabel(c,hs,'labelspacing',10000);
end
if cont(2), cc=prm.ccolor{n}; else cc='none'; end
cm=get(gcf,'colormap'); nc=size(cm,1); s=(cs(end)-cs(1))/(nc-1);
h=[]; i=1;
for hh=hs'
    m=floor((get(hh,'userdata')-cs(1))/s)+1; if m<1, m=1; elseif m>nc, m=nc; end
    if cont(1), fc=cm(m,:); else fc='none'; end
    if cont(2)&strcmp(cc,'none'), ccc=cm(m,:); else ccc=cc; end
    [lon,lat]=gmt('gridtoll',get(hh,'xdata'),get(hh,'ydata'),gprm);
    h=[h;gmt('mplot',lon,lat,ccc)];
end
[fn,fs]=gut('getfont','g'); if strcmp(cc,'none'), cc='k'; end
for hh=hl'
    pos=get(hh,'position'); str=get(hh,'string');
    [lon,lat]=gmt('gridtoll',pos(1),pos(2),gprm);
    h=[h;gmt('mtext',lon,lat,str,'color',cc,'horizontal','center','fontname',fn,...
             'fontsize',fs)];
end
delete([hs;hl])
zl=zlim;
if prm.zpos(2)<=2, gmt('setz',h,zl(prm.zpos(2))); else gmt('setz',h,zl(2)); end
if cont(1)
    caxis(cs([1,end]))
    h=[h;surf(xs,ys,zs,vs,'edgecolor','none','facecolor','interp',...
              'facelighting','phong','diffusestrength',0,'ambientstrength',1,...
              'specularstrength',0,'specularexponent',100,...
              'specularcolorreflectance',0)];
    h=[h;lightangle(0,90)];
end

% plot mesh grid --------------------------------------------------------------
function h=PlotMeshGrid(xs,ys,zs,vs,cs,prm)
i=1:prm.mint:size(xs,2); j=1:prm.mint:size(ys,1);
caxis(cs([1,end]))
h=mesh_v6(xs(j,i),ys(j,i),zs(j,i),vs(j,i),'facecolor','none');
if ~strcmp(prm.mcolor,'none'), set(h,'edgecolor',prm.mcolor); end

% plot surface ----------------------------------------------------------------
function h=PlotSurface(xs,ys,zs,vs,cs,prm)
if strcmp(prm.fcolor,'none'), c='interp'; else c=prm.fcolor; end
caxis(cs([1,end]))
h=surf(xs,ys,zs,vs,'edgecolor','none','facecolor',c,'facelighting','phong',...
       'diffusestrength',prm.lstr(1),'ambientstrength',prm.lstr(2),...
       'specularstrength',prm.fparam(1),'specularexponent',prm.fparam(2),...
       'specularcolorreflectance',prm.fparam(3));
h=[h;lightangle(prm.light(1),prm.light(2))];

% plot vector arrow -----------------------------------------------------------
function h=PlotVecArrow(xs,xe,xn,ys,ye,yn,zs,vs,u,v,cs,prm)
h=[];
scale=prm.vlen*0.001;
cm=get(gcf,'colormap'); nc=size(cm,1); s=(cs(end)-cs(1))/(nc-1);
ns=floor((vs-cs(1))/s)+1; ns(ns>nc)=nc;
for n=1:nc
    [j,i]=find(ns==n&zs>=0);
    k=find(mod(i,prm.vint)==0&mod(j,prm.vint)==0); i=i(k); j=j(k);
    x=zeros(length(j),1); y=x; uv=zeros(2,length(j));
    for m=1:length(j)
        x(m)=xs(j(m),i(m)); y(m)=ys(j(m),i(m));
        eu=[xe(j(m),i(m));ye(j(m),i(m))]; eu=eu/norm(eu);
        ev=[xn(j(m),i(m));yn(j(m),i(m))]; ev=ev/norm(ev);
        e=[eu,ev]*[u(j(m),i(m));v(j(m),i(m))];
        uv(:,m)=vs(j(m),i(m))*e/norm(e);
    end
    if ~isempty(x)
        if strcmp(prm.vcolor,'none'), c=cm(n,:); else c=prm.vcolor; end
        hh=quiver_v6(x,y,scale*uv(1,:)',scale*uv(2,:)',0); set(hh,'color',c);
        h=[h;hh];
    end
end
zl=zlim;
if prm.zpos(2)<=2, gmt('setz',h,zl(prm.zpos(2))); else gmt('setz',h,zl(2)); end

% plot stream line -------------------------------------------------------------
function h=PlotStreamLine(vec,prm,gprm)
hs=[];
[xs,ys]=meshgrid(1:prm.sint:size(vec,2),1:prm.sint:size(vec,1));
if prm.sdir==0|prm.sdir==2
    hs=[hs;streamline(vec(:,:,1),-vec(:,:,2),xs,ys,[0.1,prm.slen])];
end
if prm.sdir==1|prm.sdir==2
    hs=[hs;streamline(-vec(:,:,1),vec(:,:,2),xs,ys,[0.1,prm.slen])];
end
h=[];
for hh=hs'
    [lon,lat]=gmt('gridtoll',get(hh,'xdata'),get(hh,'ydata'),gprm);
    h=[h;gmt('mplot',lon,lat,prm.scolor)]; delete(hh);
end
zl=zlim;
if prm.zpos(2)<=2, gmt('setz',h,zl(prm.zpos(2))); else gmt('setz',h,zl(2)); end

% callback on cursor move ------------------------------------------------------
function cbBtnMove
if gcf~=gcbf, return, end
data=get(gcf,'userdata'); prm=data.prm; if data.lock|isempty(data.gprm), return, end
p=get(gca,'currentpoint');
[lon,lat]=gmt('xytoll',p(1,1),p(1,2));
[xl,yl]=gmt('getlim'); s1='';
if p(1,1)<xl(1)|xl(2)<p(1,1)|p(1,2)<yl(1)|yl(2)<p(1,2)|isnan(lon)|isnan(lat)
   gut('setinv','pos');
else
   gut('setvis','pos');
   h=gut('getchk','time'); if ~isempty(h), n=str2num(h{1}); else n=[]; end
   h=gut('getchk','height');
   if ~isempty(h)
       data1=data.data1{n}; data2=data.data2{n};
       if ~isempty(data1)
           m=find(strcmp(data.height1,h{1}(3:end)));
           s1=GetValue(prm.type{1},lat,lon,data1(:,:,m,:),data.gprm);
       end
       if ~isempty(data2)
           m=find(strcmp(data.height2,h{2}(3:end)));
           s1=[s1,' : ',GetValue(prm.type{2},lat,lon,data2(:,:,m,:),data.gprm)];
       end
   end
   we='W E'; sn='S N'; 
   s2=sprintf('%5.2f%c %6.2f%c',abs(lat),sn(sign(lat)+2),abs(lon),we(sign(lon)+2));
   gut('setstring','pos',[s2,' : ',s1]);
end

% grid value -------------------------------------------------------------------
function s=GetValue(type,lat,lon,data,gprm)
[x,y]=gmt('lltogrid',lon,lat,gprm);
if strcmp(type,'wind')
   u=interp2(1:gprm.nx,1:gprm.ny,double(data(:,:,1,1)),x,y);
   v=interp2(1:gprm.nx,1:gprm.ny,double(data(:,:,1,2)),x,y);
   z=sqrt(u^2+v^2);
   [u,v]=gmt('gridtouv',x,y,u,v,gprm);
   d=atan2(-u,-v)*180/pi; if d<0, d=d+360; end
   if ~isnan(z), s=sprintf('%.1f(%3.0f)',z,d); else s=''; end
else
   z=interp2(1:gprm.nx,1:gprm.ny,double(data),x,y);
   if ~isnan(z), s=sprintf('%.1f',z); else s=''; end
end

% read data --------------------------------------------------------------------
function ReadData
data=get(gcf,'userdata'); prm=data.prm;
[prm.td,prm.ts]=caltomjd(prm.tstart); [tn,prm.te]=caltomjd(prm.tend);
prm.te=prm.te+(tn-prm.td)*86400; if isempty(prm.te<=prm.ts), return, end
data.lock=1; set(gcf,'userdata',data)
data.time=[]; data.ft=[]; data.gprm=[]; data.data1={}; data.data2={};
data.height1={}; data.height2={};
gut('newmsgbar','msgbar','',[480,50],[0,0.7,1],1);
data.time=[]; data.ft=[]; t=prm.ts; tt=0;
while t+tt<=prm.te
    data.time=[data.time,t]; data.ft=[data.ft,tt];
    if prm.ft==0|tt+prm.ft>=prm.tint, t=t+prm.tint*3600; tt=0; else tt=tt+prm.ft; end
end
abort=0;
for n=1:length(data.time)
    msg=['reading gpv data : ',tstr(prm.td,data.time(n),data.ft(n))];
    abort=gut('setmsgbar','msgbar',msg,(n-1)/length(data.time));
    if abort, break, end
    [data1,gprm1,layer1]=...
        readgpv(prm.td,data.time(n),prm.type{1},prm.dirs.gpv,prm.src,prm.height(1),data.ft(n));
    [data1,height1]=ExtractData(prm.type{1},data1,layer1,prm.src);
    if ~isempty(gprm1)&isempty(data.gprm), data.gprm=gprm1; end
    if ~isempty(prm.type{2})
        [data2,gprm2,layer2]=...
            readgpv(prm.td,data.time(n),prm.type{2},prm.dirs.gpv,prm.src,prm.height(2),data.ft(n));
        if ~isempty(gprm2)&isempty(data.gprm), data.gprm=gprm2; end
        if ~isempty(data2)&gprm2.nx==gprm1.nx&gprm2.ny==gprm1.ny
            [data2,height2]=ExtractData(prm.type{2},data2,layer2,prm.src);
        else
            data2=[]; height2={};
        end
    else
        data2=[]; height2={};
    end
    data.data1={data.data1{:},data1};
    data.data2={data.data2{:},data2};
    data.height1={data.height1{:},height1{:}};
    data.height2={data.height2{:},height2{:}};
end
[data.height1,i]=unique(data.height1);
[data.height2,j]=unique(data.height2);
for n=1:length(data.data1), data.data1{n}=data.data1{n}(:,:,i,:); end
for n=1:length(data.data2), data.data2{n}=data.data2{n}(:,:,j,:); end
gut('closemsgbar','msgbar');
ti=[prm.dirs.gpv,' (',prm.src,')'];

if ~isempty(data.gprm)
    if data.gprm.type==3
        prm.map.proj='lambert'; prm.map.base=[data.gprm.lov,0];
    end
end
menu={}; val={};
for n=1:length(data.height1)
    menu={menu{:},data.height1{n}}; val={val{:},['1_',data.height1{n}]};
end
if ~isempty(data.data2)&~isempty(data.data2{1})
    menu={menu{:},'-'}; val={val{:},''};
    for n=1:length(data.height2)
        menu={menu{:},data.height2{n}}; val={val{:},['2_',data.height2{n}]};
    end
end
delete(gut('geth','height'));
gut('newmenu','height','&Height',menu,val,[mfilename,' cbHeight']);
if ~isempty(data.height1), gut('chkmenu','height',['1_',data.height1{end}]); end
if ~isempty(data.data2)&~isempty(data.data2{1})&~isempty(data.height2)
    gut('chkmenu','height',['2_',data.height2{end}]);
end
delete(gut('geth','time'));
menu={}; val={};
for n=1:length(data.time)
    menu={menu{:},tstr(prm.td,data.time(n),data.ft(n))}; val={val{:},num2str(n)};
end
gut('newmenu','time','&Time',menu,val,[mfilename,' cbTime']);
gut('chkmenu','time','1');
for n=1:length(data.gplot), delete(data.gplot{n}); end
data.gplot={}; for n=1:length(data.time), data.gplot{n}=[]; end
data.prm=prm; data.lock=0;
set(gcf,'name',['Grid Point Values : ',ti],'userdata',data)

% extract data -----------------------------------------------------------------
function [data,height]=ExtractData(type,data,layer,src)
height={};
for m=1:length(layer)
    if layer(m)==0, s='Surface'; else s=sprintf('%4.0fhPa',layer(m)); end
    height={height{:},s};
end

% locate gui objects -----------------------------------------------------------
function LocatePos
p=get(gcf,'position');
gut('setpos','tbtn_1',[p(3)-84,p(4)-16]);
gut('setpos','tbtn_2',[p(3)-68,p(4)-16]);
gut('setpos','tbtn_3',[p(3)-52,p(4)-16]);
gut('setpos','tbtn_4',[p(3)-36,p(4)-16]);
gut('setpos','pos',[5,p(4)-40]);

% time string ------------------------------------------------------------------
function s=tstr(td,ts,ft)
t=mjdtocal(td,ts); s=sprintf('%04d/%02d/%02d %02d:%02d T=%2d',t(1:5),ft);

% matlab 6 compatible plots ----------------------------------------------------
function [c,h]=contour_v6(varargin)
v=version; v=str2num(v(1));
if v<=6, [c,h]=contour(varargin{:}); else [c,h]=contour('v6',varargin{:}); end

function [c,h]=contourf_v6(varargin)
v=version; v=str2num(v(1));
if v<=6, [c,h]=contourf(varargin{:}); else [c,h]=contourf('v6',varargin{:}); end

function h=mesh_v6(varargin)
v=version; v=str2num(v(1));
if v<=6, h=mesh(varargin{:}); else h=mesh('v6',varargin{:}); end

function h=quiver_v6(varargin)
v=version; v=str2num(v(1));
if v<=6, h=quiver(varargin{:}); else h=quiver('v6',varargin{:}); end
