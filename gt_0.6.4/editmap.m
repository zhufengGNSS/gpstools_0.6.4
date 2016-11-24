function [prm,ok]=editmap(title,prm,opt)
%-------------------------------------------------------------------------------
% [system] : GpsTools
% [module] : map area editor dialog
% [func]   : map area editor dialog
% [argin]  : title = dialog title
%            (prm) = input map parameters
%            (opt) = edit options
%                opt(1) = 0:fix map projection
%                opt(2) = 0:fix center lon/lat
%                opt(3) = 0:fix base lon
%                opt(4) = 0:fix scale
%                opt(5) = 0:fix colors
%                opt(6) = 0:fix grid/label interval
% [argout] : prm   = output map parameters
%            ok    = exit status (1:ok,0:cancel)
% [note]   : modal dialog
% [version]: $Revision: 12 $ $Date: 2008-11-25 10:02:15 +0900 (火, 25 11 2008) $
%            Copyright(c) 2004-2006 by T.Takasu, all rights reserved
% [history]: 04/12/05  0.1  new
%            05/06/28  0.2  add map list
%            05/10/04  0.3  add mollweide projection
%-------------------------------------------------------------------------------
if nargin<1, title=''; elseif strncmp(title,'cb',2), feval(title); return, end
if nargin<2
    prm.proj='eq';
    prm.cent=[0,0];
    prm.base=[0,0];
    prm.scale=1;
    prm.gint=[10,30];
    prm.color={'w','none','b','k'};
end
if nargin<3, opt=[1,1,1,1,1,1,1]; end
[prm,ok]=EditMapArea(title,prm,opt);

% map area editor dialog -------------------------------------------------------
function [prm,ok]=EditMapArea(title,prm,opt)
strs={'Equidistant Cylindrical','Miller Cylindrical','Mercator','Orthographic',...
      'Lambert Conformal','Polar Stereo N','Polar Stereo S','Mollweide'};
vals={'eq','miller','mercat','ortho','lambert','polarn','polars','mollw'};
prm2={
'n','Center Lat',0
'n','Center Lon',0
'n','Base Lon',  0
'n','Scale',     1
'n','Grid Int.', 15
'n','Label Int.',30
};
prm3={
'c','Land Color',[0.5,0.5,0.5]
'c','Sea Color', 'w'
'c','Coast Line','b'
'c','Grid Line', 'k'
};
h=gut('newdlg','',title,[590,345]);
set(gut('newtext','pos',[377,325,90,14],'',2),'visible','off','backgroundcolor','w');
gut('newtext','',[475,325,114,20],'Map Projection',2);
gut('newpopm','prm1',[475,310,114,23],strs,vals,[mfilename,' cbUpdate']);
gut('newprms','prm2',[475,180,113,22,56],prm2);
gut('newprms','prm3',[475,92,113,22,56],prm3);
gut('newslid','scale',[2,3,15,339],[0,3],[mfilename,' cbScale']);
if ~opt(1), gut('setdis','prm1'); end
if ~opt(2), gut('setdis',{'prm2_1','prm2_2'}); end
if ~opt(3), gut('setdis','prm2_3'); end
if ~opt(4), gut('setdis',{'prm2_4','scale'}); end
if ~opt(5), gut('setdis',{'prm2_5','prm2_6'}); end
if ~opt(6), gut('setdis',{'prm3_1','prm3_2'}); end
if ~opt(7), gut('setdis',{'prm3_3','prm3_4'}); end
gut('newbtnh','',[474,68,115,22],{'Update','Reset'},...
    {[mfilename,' cbUpdate'],[mfilename,' cbReset']});
gut('newbtnv','',[475,2,113,66],{'Map List',' OK ',' Cancel '},...
    {[mfilename,' cbList'],[mfilename,' cbOK'],[mfilename,' cbCancel']});
prm.p0=[]; DrawMap(prm); SetPrm(prm);
set(h,'doublebuffer','on','windowbuttondownfcn',[mfilename,' cbBtnDown'],...
    'windowbuttonmotionfcn',[mfilename,' cbBtnMove'],...
    'windowbuttonupfcn',[mfilename,' cbBtnUp']);
ok=gut('waitok');
if ok, prm=GetPrm; close; end
rmfield(prm,'p0');

% set map parameters -----------------------------------------------------------
function SetPrm(prm)
if isempty(prm), return, end
prm.cent=round(prm.cent*1E3)/1E3;
prm.scale=round(prm.scale*1E3)/1E3;
gut('setsel','prm1',prm.proj);
gut('setprms','prm2',prm.cent(2),prm.cent(1),prm.base(1),prm.scale,prm.gint(1),...
    prm.gint(2));
gut('setprms','prm3',prm.color{1},prm.color{2},prm.color{3},prm.color{4});
gut('setval','scale',log10(prm.scale));
set(gca,'userdata',prm);

% get map parameters -----------------------------------------------------------
function prm=GetPrm
prm=get(gca,'userdata');
prm.proj=gut('getsel','prm1');
[prm.cent(2),prm.cent(1),prm.base(1),prm.scale,prm.gint(1),prm.gint(2)]=...
    gut('getprms','prm2');
[prm.color{1},prm.color{2},prm.color{3},prm.color{4}]=gut('getprms','prm3');
prm.cent(1)=gmt('normlon',prm.cent(1));
prm.cent(2)=min(max(prm.cent(2),-90),90);
prm.base(1)=gmt('normlon',prm.base(1));
prm.scale=min(max(prm.scale,1),1000);

% draw map ---------------------------------------------------------------------
function DrawMap(prm)
[fn,fs]=gut('getfont','g');
cla, axis off
gmt('mmap','proj',prm.proj,'cent',prm.cent,'base',prm.base,'scale',prm.scale,...
    'pos',[0.024,0.008,0.78,0.98],'color','none','fontname',fn,'fontsize',fs)
gmt('mcoast','lcolor',prm.color{1},'scolor',prm.color{2},'ccolor',prm.color{3});
gmt('mgrid','gint',prm.gint(1),'lint',prm.gint(2),'color',prm.color{4});
drawnow

% show current position --------------------------------------------------------
function ShowPos(p)
[lon,lat]=gmt('xytoll',p(1),p(2)); [xl,yl]=gmt('getlim');
if p(1,1)<xl(1)|xl(2)<p(1)|p(2)<yl(1)|yl(2)<p(2)|isnan(lon)|isnan(lat)
    gut('setinv','pos');
else
    we='W E'; sn='S N'; 
    s=sprintf('%5.2f%c %6.2f%c',abs(lat),sn(sign(lat)+2),abs(lon),we(sign(lon)+2));
    gut('setstring','pos',s); gut('setvis','pos');
end

% callback on ok/cancel --------------------------------------------------------
function cbOK,     set(gcf,'userdata',1);
function cbCancel, set(gcf,'userdata',0);

% callback on map list ---------------------------------------------------------
function cbList
list={}; prms={};
[path,f]=fileparts(which(mfilename));
file=fullfile(path,'settings','prm_editmap.mat');
if exist(file), load(file); end
[list,prms,sel]=editprms('Map List',list,prms,'',GetPrm);
save(file,'list','prms')
if isempty(prms)|isempty(sel), return, end
p=prms{sel}; prm=GetPrm; prm.proj=p.proj; prm.cent=p.cent; prm.base=p.base;
prm.scale=p.scale; prm.gint=p.gint;
DrawMap(prm); SetPrm(prm);

% callback on update parameters ------------------------------------------------
function cbUpdate
prm=GetPrm; DrawMap(prm); SetPrm(prm);

% callback on change scale -----------------------------------------------------
function cbScale
prm=GetPrm; prm.scale=10^gut('getval','scale'); SetPrm(prm); DrawMap(prm);

% callback on reset ------------------------------------------------------------
function cbReset
prm=GetPrm; prm.cent=[0,0]; prm.scale=1; SetPrm(prm); DrawMap(prm);

% callback on button down ------------------------------------------------------
function cbBtnDown
prm=get(gca,'userdata'); p=get(gca,'currentpoint'); 
[lon,lat]=gmt('xytoll',p(1,1),p(1,2)); if isnan(lon)|isnan(lat), return, end
switch get(gcf,'selectiontype')
case 'normal'
    prm.p0=p(1,1:2); set(gca,'userdata',prm);
    set(gcf,'pointer','fleur');
case 'open'
    prm.cent=[lon,lat]; SetPrm(prm); DrawMap(prm); % double click
end

% callback on button move ------------------------------------------------------
function cbBtnMove
prm=get(gca,'userdata'); p=get(gca,'currentpoint'); 
if ~isempty(prm.p0)
    cent=MapCent(p(1,1:2)-prm.p0,prm);
    if all(~isnan(cent)), prm.cent=cent; DrawMap(prm); end
else
    ShowPos(p(1,:))
end

% callback on button up -------------------------------------------------------
function cbBtnUp
prm=get(gca,'userdata'); p=get(gca,'currentpoint'); 
if strcmp(get(gcf,'selectiontype'),'normal')&~isempty(prm.p0)
    cent=MapCent(p(1,1:2)-prm.p0,prm); if all(~isnan(cent)), prm.cent=cent; end
    prm.p0=[]; SetPrm(prm); DrawMap(prm); 
    set(gcf,'pointer','arrow');
end

% get map center parameter -----------------------------------------------------
function cent=MapCent(dp,prm)
switch prm.proj
case 'ortho', cent=prm.cent-dp*90/prm.scale;
otherwise, [lon,lat]=gmt('xytoll',-dp(1),-dp(2),prm); cent=[lon,lat];
end
cent(1)=gmt('normlon',cent(1));
cent(2)=min(max(cent(2),-90),90);
