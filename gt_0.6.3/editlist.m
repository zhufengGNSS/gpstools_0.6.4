function [list,ok,label]=editlist(title,list,varargin)
%-------------------------------------------------------------------------------
% [system] : GpsTools
% [module] : satellite/station list editor dialog
% [func]   : satellite/station list editor dialog
% [argin]  : title = dialog title
%            list  = list
%            opts  = option
%                    'rcv'    : station list
%                    'noedit' : disable edit
% [argout] : list  = list
%            ok    = exit status (1:ok,0:cancel)
%            label = label strings
% [note]   : modal dialog
% [version]: $Revision: 4 $ $Date: 06/07/20 13:30 $
%            Copyright(c) 2004-2006 by T.Takasu, all rights reserved
% [history]: 04/12/05  0.1  new
%            05/06/28  0.2  add station map view
%            06/07/06  0.3  delete argin prm
%-------------------------------------------------------------------------------
if strncmp(title,'cb',2), feval(title)
else [list,ok,label]=EditList(title,list,varargin{:}); end

% edit list body ---------------------------------------------------------------
function [list,ok,label]=EditList(title,list,varargin)
label={}; selrcv=0; noedit=0; n=1;
while n<=nargin-2
    switch varargin{n}
    case 'rcv',    selrcv=1; n=n+1;
    case 'noedit', noedit=1; n=n+1;
    otherwise, n=n+1; end
end
h=gut('newdlg','edlg',title,[535,400]); set(h,'renderer','painters');
if selrcv
    prm=prm_gpsrcvs; grp=unique(prm(:,2));
else
    prm=prm_gpssats; grp={''};
end
set(h,'userdata',prm);
gut('newtext','text1',[10,375,230,20],'',2);
gut('newpopm','grp', [301,376,228,23],grp,{},[mfilename,' cbChgGrp']);
if selrcv, gut('setsel','grp',prm{1,2}); end
gut('newlist','list',[10,30,230,348],{},{},2,[mfilename,' cbListClk']);
gut('newlist','prms',[300,30,230,348],{},{},2,[mfilename,' cbPrmClk']);
gut('newbtnv','abtn',[245,168,50,72],{'<<','>>','Clear'},...
    {[mfilename,' cbAdd'],[mfilename,' cbDel'],[mfilename,' cbClear']});
gut('newbtn','ebtn',[245,90,50,23],'Edit...',[mfilename,' cbEdit']);
gut('newbtnh','btns',[10,4,225,23],{'Load...','Save...','List...','Input...'},...
    {[mfilename,' cbLoad'],[mfilename,' cbSave'],[mfilename,' cbPrms'],...
     [mfilename,' cbInput']});
gut('newbtn','mapb',[240,4,60,23],'MapView',[mfilename,' cbMap']);
gut('newslid','scale',[2,31,14,344],[0,3],[mfilename,' cbScale']);
gut('setudata','btns_1',selrcv);
gut('newokcancelbtn','',[330,4,180,23]);
if noedit
    gut('setdis',{'prms','grp','add','del','btns_1','btns_2','btns_3','btns_4'});
end
if ~selrcv, gut('setdis',{'mapb','scale'}); end
cbChgGrp
SetList(list,selrcv);
prm.proj='miller'; prm.cent=[0,0]; prm.scale=1; prm.gint=[15,30]; prm.p0=[];
gut('setudata','mapb',prm); gut('setinv','scale');
set(h,'doublebuffer','on','windowbuttondownfcn',[mfilename,' cbBtnDown']);
UpdateLabel;
ok=gut('waitok'); if ~ok, return, end
list=gut('getudata','list');
label=gut('getstring','list');
close

% callback on load -------------------------------------------------------------
function cbLoad
selrcv=get(gcbo,'userdata');
patt={'*.txt','text (*.txt)';'prm_*.mat','settings (prm_*.mat)';'*.*','all (*.*)'};
[file,path]=uigetfile(patt,'Load List');
if file==0, return, end
[d,f,ext]=fileparts(file);
switch ext
case '.mat'
    list={}; prm=[];
    load([path,file])
    if ~isempty(prm)
        if selrcv, list=prm.rcvs; else list=prm.sats; end
    else
        for n=1:length(list), list{n}=strtok(list{n}); end
    end
case '.m'
    wd=pwd; cd(path); prm=feval(f); cd(wd);
    if selrcv, list=prm.rcvs; else list=prm.sats; end
case '.txt'
    list=textread([path,file],'%s%*[^\n]');
    list=sort(upper(list));
end
SetList(list,selrcv);
UpdateLabel;
if strcmp(gut('getstring','mapb'),'ListView'), DrawMap; end

% callback on save -------------------------------------------------------------
function cbSave
selrcv=get(gcbo,'userdata');
patt={'*.txt','text (*.txt)';'*.*','all (*.*)'};
[file,path]=uiputfile(patt,'Save List');
if file==0, return, end
[d,file,ext]=fileparts(file);
f=fopen([path,file,ext],'wt'); if f<0, return, end
[str,val]=gut('getlist','list');
for n=1:length(val), fprintf(f,'%s\n',val{n}); end
fclose(f);

% callback on parameters list --------------------------------------------------
function cbPrms
list={}; prms={}; selrcv=gut('getudata','btns_1');
if selrcv, s='Receiver'; f='rcvs'; else s='Satellite'; f='sats'; end
[dirs,fs]=fileparts(which(mfilename));
file=fullfile(dirs,'settings',['prm_editlist_',f,'.mat']);
if exist(file), load(file); end
[label,slist]=gut('getlist','list');
[list,prms,sel]=editprms([s,' List'],list,prms,'',slist);
save(file,'list','prms')
if ~isempty(prms)&~isempty(sel)
    SetList(prms{sel},selrcv);
    if strcmp(gut('getstring','mapb'),'ListView'), DrawMap; end
end

% callback on parameters input -------------------------------------------------
function cbInput
gut('newdlg','','Input Item',[250,55]);
h=gut('newedit','inp',[5,30,240,20],'',1);
set(h,'callback',[mfilename,' cbInputOk']);
gut('newokcancelbtn','',[65,4,180,23]);
if ~gut('waitok'), return, end
inp=gut('getstring','inp');
close
if ~isempty(inp), AddList({inp},{inp}); end

function cbInputOk, set(gcf,'userdata',1)

% callback on map --------------------------------------------------------------
function cbMap
if strcmp(get(gcbo,'string'),'MapView')
    gut('setinv',{'prms','list','abtn_1','abtn_2','abtn_3','ebtn'});
    gut('setvis','scale');
    DrawMap;
    set(gcbo,'string','ListView');
else
    gut('setvis',{'prms','list','abtn_1','abtn_2','abtn_3','ebtn'});
    gut('setinv','scale');
    cla;
    set(gcbo,'string','MapView');
end

% callback on edit -------------------------------------------------------------
function cbEdit
if strcmp(get(gut('geth','mapb'),'enable'),'on'), prmfile='rcvs_params.txt';
else prmfile='sats_params.txt'; end
[dirs,f]=fileparts(which(mfilename));
file=fullfile(dirs,'data',prmfile);
if exist(file)
    p=gut('getgutprm');
    [s,w]=dos(['"',p.editor,'" "',file,'" &;exit']);
    if s, gut('errdlg',w); end
end

% callback on button down ------------------------------------------------------
function cbBtnDown
if strcmp(gut('getstring','mapb'),'MapView'), return, end
prm=gut('getudata','mapb');
p=get(gca,'currentpoint'); 
[lon1,lat1]=gmt('xytoll',p(1,1),p(1,2)); if isnan(lon1)|isnan(lat1), return, end
if strcmp(get(gcf,'selectiontype'),'open') % double click
    if -90<lat1&lat1<90
        prm.cent=[lon1,lat1]; gut('setudata','mapb',prm); DrawMap;
    end
else
    rect=rbbox; p=get(gca,'currentpoint'); 
    [lon2,lat2]=gmt('xytoll',p(1,1),p(1,2));
    if isnan(lon2)|isnan(lat2)|(lon1==lon2&lat1==lat2), return, end
    add=strcmp(get(gcf,'selectiontype'),'alt');
    SelectRcvs([lon1,lon2],[lat1,lat2],prm.cent(1),add);
end

% callback on change scale -----------------------------------------------------
function cbScale
prm=gut('getudata','mapb');
prm.scale=10^gut('getval','scale');
gut('setudata','mapb',prm);
DrawMap;

% draw map ---------------------------------------------------------------------
function DrawMap
prm=gut('getudata','mapb');
[fn,fs]=gut('getfont','g');
cla, axis off
gmt('mmap','proj',prm.proj,'cent',prm.cent,'base',[0,0],'scale',prm.scale,...
    'pos',[0.03,0.06,0.965,0.89],'color','none','fontname',fn,'fontsize',fs)
gmt('mcoast','lcolor','w','scolor','w','ccolor',[0.7,0.7,0.7]);
gmt('mgrid','gint',prm.gint(1),'lint',prm.gint(2),'color',[0.5,0.5,0.5]);
[str,rcvs]=gut('getlist','prms'); DrawRcvs(rcvs,[0.7,0.7,0.7],0)
[str,rcvs]=gut('getlist','list'); DrawRcvs(rcvs,'r',0)
gut('setval','scale',log10(prm.scale));
drawnow

% draw station pos -------------------------------------------------------------
function DrawRcvs(rcvs,color,opt)
if isempty(rcvs), return, end
[fn,fs]=gut('getfont','g');
xl=get(gca,'xlim'); yl=get(gca,'ylim');
posr=readpos(0,0,rcvs,'','approx');
for n=1:length(rcvs), gpos(n,:)=eceftogeod(posr(1,:,n)'); end
[x,y,z]=gmt('lltoxy',gpos(:,2),gpos(:,1));
i=find(xl(1)<=x&x<=xl(2)&yl(1)<=y&y<=yl(2)&z>=0);
plot(x(i),y(i),'.','color',color);
if opt
    for n=i'
         gmt('mtext',gpos(n,2),gpos(n,1),rcvs{n},'horizontal','center',...
             'vertical','top','fontname',fn,'fontsize',fs);
    end
end

% select station select --------------------------------------------------------
function SelectRcvs(lons,lats,lonc,add)
[label,rcvs]=gut('getlist','prms');
posr=readpos(0,0,rcvs,'','approx'); lon=zeros(length(rcvs),1); lat=lon;
for n=1:length(rcvs)
    p=eceftogeod(posr(1,:,n)'); lon(n)=gmt('normlon',p(2)-lonc); lat(n)=p(1);
end
lons=sort(gmt('normlon',lons-lonc)); lats=sort(lats);
i=find(lons(1)<=lon&lon<=lons(2)&lats(1)<=lat&lat<=lats(2));
if add
    AddList(rcvs(i),label(i));
else
    gut('setlist','list',label(i),rcvs(i));
end
DrawMap; UpdateLabel;

% callback on group change -----------------------------------------------------
function cbChgGrp
prm=get(gcf,'userdata');
grp=gut('getsel','grp');
labels={}; rcvs={}; m=0;
for n=1:size(prm,1)
    if isempty(grp), m=m+1; labels{m}=prm{n,1}; rcvs{m}=prm{n,1};
    elseif strcmp(prm{n,2},grp)
        m=m+1; labels{m}=[prm{n,1},' ',prm{n,3}]; rcvs{m}=prm{n,1};
    end
end
gut('setlist','prms',labels,rcvs);
mprm=gut('getudata','mapb');
if strcmp(grp,'GSI')
    mprm.cent=[135,35]; mprm.scale=8.5; mprm.gint=[5,10];
else
    mprm.cent=[0,0]; mprm.scale=1; mprm.gint=[15,30];
end
gut('setudata','mapb',mprm);
if ~strcmp(gut('getstring','mapb'),'MapView'), DrawMap; end

% callback on add select list --------------------------------------------------
function cbAdd
[rcvs,label]=gut('getsel','prms');
AddList(rcvs,label);

% callback on delete select list -----------------------------------------------
function cbDel
DelList(gut('getsel','list'));

% callback on clear list -------------------------------------------------------
function cbClear
gut('setlist','list',{},{});
UpdateLabel;

% callback on Parameter Clik ---------------------------------------------------
function cbPrmClk
if ~strcmp(get(gcf,'selectiontype'),'open'), return, end
[rcvs,label]=gut('getsel','prms');
AddList(rcvs,label);

% callback on List Clik --------------------------------------------------------
function cbListClk
if ~strcmp(get(gcf,'selectiontype'),'open'), return, end
[rcvs,label]=gut('getsel','list');
DelList(rcvs);

% add list ---------------------------------------------------------------------
function AddList(rcvs,label)
[str,val]=gut('getlist','list');
str={str{:},label{:}};
val={val{:},rcvs{:}};
[val,i]=unique(val); str=str(i);
gut('setlist','list',str,val);
gut('setsel','list',rcvs);
UpdateLabel;

% delete list ------------------------------------------------------------------
function DelList(rcvs)
[str,val]=gut('getlist','list'); i=[];
for n=1:length(rcvs), i=[i,find(strcmp(val,rcvs{n}))]; end
str(i)=[]; val(i)=[];
gut('setlist','list',str,val);
UpdateLabel;

% set list ---------------------------------------------------------------------
function SetList(list,selrcv)
prm=get(gcf,'userdata');
labels={}; vals={}; m=0;
for n=1:length(list)
    i=strcmp(prm(:,1),list{n});
    if selrcv&~isempty(i), label=[list{n},' ',prm{i,3}]; else label=list{n}; end
    m=m+1;
    labels{m}=label;
    vals{m}=list{n};
end
gut('setlist','list',labels,vals);
UpdateLabel;

% update label -----------------------------------------------------------------
function UpdateLabel
str=gut('getlist','list');
gut('setstring','text1',sprintf('Selected %s (%d)',get(gcf,'name'),length(str)));
