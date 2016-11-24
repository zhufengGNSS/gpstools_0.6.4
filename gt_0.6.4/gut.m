function varargout=gut(varargin)
%-------------------------------------------------------------------------------
% [system] : GpsTools
% [module] : gui common libraries for GpsTools
% [func]   : gui common libraries for GpsTools
% [argin]  : function,args
%              'initgut'(,prm)               : initialize gut
%              'setfonts',type,fname,fsize   : set font name/size
%              'getfonts',type               : get font name/size
%              'newfig',tag,title,siz(,cb)   : generate figure
%              'newfigg',tag,title,siz       : generate graph figure
%              'showfig',tag                 : show figure
%              'hidefig',tag                 : hide figure
%              'newdlg',tag,title,siz(,pos)  : generate modal dialog
%              'newtext',tag,pos,str(,align) : generate text
%              'newtextv',tag,pos,strs(,align) : generate texts column
%              'newtexth',tag,pos,strs(,align) : generate texts row
%              'newedit',tag,pos,str(,align) : generate text edit
%              'newfrm',tag,pos              : generate frame
%              'newslid',tag,pos,ext,cbs     : generate slider
%              'newchk',tag,pos,lable(,cb)   : generate check
%              'geth',tag                    : get object handle
%              'existh',h                    : check handle existing
%              'getval',tag                  : get value
%              'setval',tag,val              : set value
%              'getstring',tag               : get string
%              'setstring',tag,str           : set string
%              'getnum',tag                  : get numerical value
%              'setnum',tag,num              : set numerical value
%              'getudata',tag                : get user data
%              'setudata',tag,udata          : set user data
%              'setena',tags                 : set enable
%              'setdis',tags                 : set disable
%              'setenah',h                   : set enable by handle
%              'setdish',h                   : set disable by handle
%              'setenaall',h                 : set enable all
%              'setdisall',h                 : set disable all
%              'setvis',tags                 : set visible
%              'setinv',tags                 : set invisible
%              'getpos',tag                  : get position
%              'setpos',tag,pos              : set position
%              'getsize',tag                 : get size
%              'setsize',tag,size            : set size
%              'newbtn',tag,pos,label,cb     : generate button
%              'newbtnv',tag,pos,labels,cbs  : generate buttons row
%              'newbtnh',tag,pos,labels,cbs  : generate buttons column
%              'newtoggle',tag,pos,cb        : generate toggle
%              'newtable',tag,pos,n,m,value  : generate table
%              'gettable',tag                : get table values
%              'settable',tag,value          : set table values
%              'newmenu',tag,top,labels,cbs  : generate menu
%              'newcmenu',tag,pos,labels,cbs : generate context menu
%              'showmenu',tag                : show menu
%              'hidemenu',tag                : hide menu
%              'chkmenu',tag(,val)           : check menu
%              'unchkmenu',tag(,val)         : uncheck menu
%              'enamenu',h(,val)             : enable menu
%              'dismenu',h(,val)             : disable menu
%              'getchk',tag                  : get checked menu values
%              'togglechk',h                 : toggle menu checked
%              'newokbtn',tag,pos            : generate ok button
%              'newcancelbtn',tag,pos        : generate cancel button
%              'newokcancelbtn',tag,pos      : generate ok and cancel button
%              'newlist',tag,pos,labels(,vals,smode,cb) : generate list
%              'newlistf',tag,pos,labels(,vals,smode,cb) : generate fixed width font list
%              'newpopm',tag,pos,labels(,vals,cb) : generate combo box
%              'setlist',tag,labels(,vals)   : set list or combo box
%              'getlist',tag                 : get list or combo box
%              'setsel',tag,sel              : set select in list or combo box
%              'getsel',tag                  : get select in list or combo box
%              'newprog',tag,pos(,color)     : generate progress bar
%              'setprog',tag,prog            : set progress bar
%              'setprogh',h,prog             : set progress bar by handle
%              'newymd',tag,pos,ymd(,cb)     : generate date input
%              'getymd',tag                  : get date input
%              'setymd',tag,ymd              : set date input
%              'newtime',tag,pos,time(,opt)  : generate date/time input
%              'gettime',tag                 : get date/time input
%              'settime',tag,time            : set date/time input
%              'newprms',tag,pos,prms        : generate parameters input
%              'setprms',tag,prm1,prm2,...   : set parameters input
%              'getprms',tag                 : get parameters input
%              'updatefilelist',tag,dirs(,patt,order,ud,opt) : update file list
%              'getfilelist',tag,dirs        : get select in file list
%              'updatedirlist',tag,dirs      : update directory list
%              'getdirlist',tag              : get select in directory list
%              'newdedit',tag,pos,str        : generate directory input edit
%              'newfedit',tag,pos,str(,opt)  : generate file input edit
%              'newviewer',tag,title,siz,opt : generate text viewer
%              'newmsgbar',tag,title,siz(,color) : generate message and progress bar
%              'setmsgbar',tag,msg,prog      : set message and progress bar
%              'dirdlg',dirs                 : generate directory select dialog
%              'newdirdlg'                   : generate new directory dialog
%              'aboutdlg',title,siz,msgs     : generate about dialog
%              'seldlg',title,siz,labels(,vals,sel,smode) : generate select list dialog
%              'errdlg',msg                  : generate error message dialog
%              'confdlg',msg(,btns)          : generate confirmation dialog
%              'waitok'                      : wait for pushing ok/cancel button
% [argout] : argouts
% [note]   :
% [version]: $Revision: 16 $ $Date: 2008-12-12 15:49:30 +0900 (金, 12 12 2008) $
%            Copyright(c) 2004-2008 by T.Takasu, all rights reserved
% [history]: 04/12/06  0.1  new
%            05/02/23  0.2  set figure background color same as uicontrol
%            05/04/25  0.3  text viewer supports compact/compressed files
%            05/12/02  0.4  set figure renderermode to manual
%            06/03/07  0.5  set figure renderermode to auto
%            06/06/23  0.6  modify initialization of global parameters
%            08/12/05  0.7  fix bug on textviewer button placement if error (gt_0.6.4)
%                           fix bug on edit-text input on dialog
%                           fix bug on uigetfile in matlab 7.4
%                           add list font to global p_gut
%                           add option to newtime
%                           support file date format change in matlab R2007b
%                           suppress warnings
%-------------------------------------------------------------------------------
global p_gut, if isempty(p_gut), initgut, end
if nargout==0, feval(varargin{:});
else [varargout{1:nargout}]=feval(varargin{:}); end

% initialise global parameters -------------------------------------------------
function initgut(prm)
global p_gut
[path,file]=fileparts(which(mfilename));
prmfile=fullfile(path,'settings','prm_gut.mat');
if nargin>=1
    save(prmfile,'prm'); p_gut=prm; 
elseif exist(prmfile)==2
    load(prmfile); p_gut=prm;
else
    switch lower(computer)
    case 'linux'
        p_gut.editor='vi';              % text editor
        setfont('u','helvetica',8)      % ui font name and size
        setfont('b','helvetica',10)     % button font name and size
        setfont('f','courier',9)        % fixed width font name and size
        setfont('g','helvetica',8)      % graph font name and size
        setfont('l','helvetica',8)      % list font name and size
    otherwise
        p_gut.editor='notepad.exe';     % text editor
        setfont('u','Arial',8)          % ui font name and size
        setfont('b','Arial',9)          % button font name and size
        setfont('f','Courier New',9)    % fixed width font name and size
        setfont('g','Times New Roman',8) % graph font name and size
        setfont('l','MS UI Gothic',9)   % list font name and size
    end
    p_gut.viewbc='w';                   % background color of text viewer
    p_gut.viewfc='k';                   % foreground color of text viewer
    [fname,fsize]=getfont('f');
    font.FontName=fname;
    font.FontSize=fsize;
    font.FontUnits='points';
    font.FontWeight='normal';
    font.FontAngle='normal';
    p_gut.viewfont=font;                % font of text viewer
    p_gut.viewln=0;                     % show line number in text viewer
    p_gut.maxtext=20000;                % max lines of text viewer
end
p=get(0,'screensize');                  % position of text viewer
p_gut.viewpos=[(p(3)-600)/2,(p(4)-400)/2,600,400];

function prm=getgutprm
global p_gut
prm=p_gut;

function setfont(type,fname,fsize)
global p_gut
switch type
case 'u', p_gut.ufont=fname; p_gut.usize=fsize;
case 'b', p_gut.bfont=fname; p_gut.bsize=fsize;
case 'f', p_gut.ffont=fname; p_gut.fsize=fsize;
case 'g', p_gut.gfont=fname; p_gut.gsize=fsize;
case 'l', p_gut.lfont=fname; p_gut.lsize=fsize;
end

function [fname,fsize]=getfont(type)
global p_gut
switch type
case 'u', fname=p_gut.ufont; fsize=p_gut.usize;
case 'b', fname=p_gut.bfont; fsize=p_gut.bsize;
case 'f', fname=p_gut.ffont; fsize=p_gut.fsize;
case 'g', fname=p_gut.gfont; fsize=p_gut.gsize;
case 'l', fname=p_gut.ufont; fsize=p_gut.usize;
end

% figure -----------------------------------------------------------------------
function h=newfig(tag,title,siz,cb,varargin)
if nargin<4, cb=''; end, if isempty(cb), cb='closereq'; end
h=NewF(tag,title,siz,'menubar','none','resize','off','closerequestfcn',cb,...
       varargin{:});
if ~isempty(h), figure(h); set(gca,'visible','off'); end

function h=newfigg(tag,title,siz,cb,varargin)
if nargin<4, cb=''; end, if isempty(cb), cb='closereq'; end
h=NewF(tag,title,siz,'color','w','doublebuffer','on','closerequestfcn',cb,...
       'renderer','painters','renderermode','manual',varargin{:});
if ~isempty(h), figure(h); set(gca,'visible','off'); end

function h=NewF(tag,title,siz,varargin)
if ~isempty(tag)
    h=findobj(0,'tag',tag);
    if ~isempty(h), set(h,'visible','on'), figure(h), h=[]; return, end
end
if isempty(gcbf), p=get(0,'screensize')+[0,-40,0,0]; else p=get(gcbf,'position'); end
if length(siz)>=4, off=siz(3:4); else off=[0,0]; end
h=figure('tag',tag,'name',title,'toolbar','none','numbertitle','off',...
         'position',[p(1)+(p(3)-siz(1))/2+off(1),p(2)+p(4)-siz(2)+off(2),siz(1:2)],...
         'paperpositionmode','auto','renderermode','auto',...
         'color',get(0,'defaultuicontrolbackgroundcolor'),varargin{:});

function showfig(tag), set(findobj(0,'tag',tag),'visible','on');
function hidefig(tag), set(findobj(0,'tag',tag),'visible','off');

% modal dialog -----------------------------------------------------------------
function h=newdlg(tag,title,siz,pos,varargin)
if nargin<4, pos=1; end
f=get(0,'currentfigure');
if ~isempty(f)&pos, p=get(f,'position'); else p=get(0,'screensize'); end
if length(siz)>=4, off=siz(3:4); else off=[0,0]; end
p=[p(1)+(p(3)-siz(1))/2+off(1),p(2)+(p(4)-siz(2))/2+off(2),siz(1:2)];
h=dialog('tag',tag,'name',title,'position',p,'handlevisibility','on',varargin{:});

% ui control -------------------------------------------------------------------
function h=newui(style,tag,pos,varargin)
global p_gut
if strcmp(style,'listbox'), fn=p_gut.lfont; fs=p_gut.lsize;
else fn=p_gut.ufont; fs=p_gut.usize; end
h=uicontrol('style',style,'tag',tag,'position',pos,'fontname',fn,'fontsize',fs,...
            'handlevisibility','off',varargin{:});

% text -------------------------------------------------------------------------
function h=newtext(tag,pos,str,align)
if nargin<4, align=1; end, a={'left','center','right'};
h=newui('text',tag,pos,'horizontal',a{align},'string',str,'backgroundcolor',...
        get(gcf,'color'));

function h=newtextv(tag,pos,strs,align)
if nargin<4, align=1; end
for n=1:length(strs)
    p=[pos(1),pos(2)+pos(4)*(length(strs)-n)-8,pos(3),23];
    h(n)=newtext([tag,'_',num2str(n)],p,strs{n},align);
end

function h=newtexth(tag,pos,strs,align)
if nargin<4, align=1; end
for n=1:length(strs)
    p=[pos(1)+pos(3)/length(strs)*(n-1),pos(2),pos(3)/length(strs),pos(4)];
    h(n)=newtext([tag,'_',num2str(n)],p,strs{n},align);
end

% edit -------------------------------------------------------------------------
function h=newedit(tag,pos,str,align)
if nargin<4, align=3; end, a={'left','center','right'};
h=newui('edit',tag,pos,'horizontal',a{align},'string',str);
set(h,'backgroundcolor','w');

% frame ------------------------------------------------------------------------
function h=newfrm(tag,pos)
h=newui('frame',tag,pos);

% slider -----------------------------------------------------------------------
function h=newslid(tag,pos,ext,cb)
h=newui('slider',tag,pos,'min',ext(1),'max',ext(2),'callback',cb,'backgroundcolor',...
        get(gcf,'color'));

% check box --------------------------------------------------------------------
function h=newchk(tag,pos,label,cb)
h=newui('checkbox',tag,pos,'string',label,'backgroundcolor',get(gcf,'color'));
if nargin>=4, set(h,'callback',cb); end

% find hiden object ------------------------------------------------------------
function h=findhobj(varargin)
s=get(0,'showhiddenhandles');
set(0,'showhiddenhandles','on');
h=findobj(varargin{:});
set(0,'showhiddenhandles',s);

% get handle -------------------------------------------------------------------
function h=geth(tag)
h=findhobj(gcf,'tag',tag);                         % obj in current figure
if isempty(h), h=findhobj(gcbf,'tag',tag); end     % obj in callback figure
if isempty(h), h=findhobj(0,'flat','tag',tag); end % figure

% check handle existing --------------------------------------------------------
function stat=existh(h)
s=get(0,'showhiddenhandles');
set(0,'showhiddenhandles','on');
stat=any(findobj==h);
set(0,'showhiddenhandles',s);

% get/set value ----------------------------------------------------------------
function setval(tag,val), set(geth(tag),'value',val);
function val=getval(tag), val=get(geth(tag),'value');

% get/set string ---------------------------------------------------------------
function str=getstring(tag), str=get(geth(tag),'string');
function setstring(tag,str), set(geth(tag),'string',str);

% get/set numerical value ------------------------------------------------------
function num=getnum(tag), num=StrToNum(getstring(tag));
function setnum(tag,num), setstring(tag,NumToStr(num));

% get/set font properties ------------------------------------------------------
function font=getfonth(h)
font.FontName=get(h,'fontname');
font.FontSize=get(h,'fontsize');
font.FontUnits=get(h,'fontunits');
font.FontWeight=get(h,'fontweight');
font.FontAngle=get(h,'fontangle');

function setfonth(h,font)
if isstruct(font), set(h,font), end

% get/set userdata -------------------------------------------------------------
function udata=getudata(tag), udata=get(geth(tag),'userdata');
function setudata(tag,udata), set(geth(tag),'userdata',udata);

% set enable/disable -----------------------------------------------------------
function h=setena(tags)
if ~iscell(tags), tags={tags}; end
h=[]; for n=1:length(tags), h=[h,geth(tags{n})]; end
setenah(h)

function h=setdis(tags)
if ~iscell(tags), tags={tags}; end
h=[]; for n=1:length(tags), h=[h,geth(tags{n})]; end
setdish(h)

function setenah(h)
for n=1:length(h), set(h(n),'enable','on'), end
updatenow

function setdish(h)
for n=1:length(h), set(h(n),'enable','off'), end
updatenow

function setenaall(h)
for hh=findhobj(h,'type','uicontrol')', set(hh,'enable','on'), end
updatenow

function setdisall(h)
for hh=findhobj(h,'type','uicontrol')', set(hh,'enable','off'), end
updatenow

% set visible/invisible --------------------------------------------------------
function h=setvis(tags)
if ~iscell(tags), tags={tags}; end
for n=1:length(tags), set(geth(tags{n}),'visible','on'), end

function h=setinv(tags)
if ~iscell(tags), tags={tags}; end
for n=1:length(tags), set(geth(tags{n}),'visible','off'), end

% get/set position -------------------------------------------------------------
function pos=getpos(tag)
p=get(geth(tag),'position'); pos=p(1:2);

function setpos(tag,pos)
h=geth(tag); p=get(h,'position'); set(h,'position',[pos,p(3:4)]); % unchange size

% get/set size -----------------------------------------------------------------
function size=getsize(tag)
p=get(geth(tag),'position'); size=p(3:4);

function setsize(tag,size)
h=geth(tag); p=get(h,'position'); set(h,'position',[p(1:2),size]); % unchange pos

% buttons ----------------------------------------------------------------------
function h=newbtn(tag,pos,label,cb)
global p_gut
if isempty(cb), cb='close'; end
h=newui('pushbutton',tag,pos,'string',label,'callback',cb,'backgroundcolor',...
        get(gcf,'color'));
set(h,'fontname',p_gut.bfont,'fontsize',p_gut.bsize);

function h=newbtnh(tag,pos,labels,cbs)
w=pos(3)/length(labels);
for n=1:length(labels)
    if iscell(cbs), cb=cbs{n}; else cb=cbs; end
    h(n)=newbtn([tag,'_',num2str(n)],[pos(1)+(n-1)*w+1,pos(2),w-1,pos(4)],...
                labels{n},cb);
end

function h=newbtnv(tag,pos,labels,cbs)
if ~iscell(labels), labels={labels}; end
hh=pos(4)/length(labels);
for n=1:length(labels)
    if iscell(cbs), cb=cbs{n}; else cb=cbs; end
    h(n)=newbtn([tag,'_',num2str(n)],[pos(1),pos(2)+pos(4)-n*hh+1,pos(3),hh-2],...
                labels{n},cb);
end

% toggle -----------------------------------------------------------------------
function h=newtoggle(tag,pos,cb)
global p_gut
if isempty(cb), cb='close'; end
h=newui('toggle',tag,pos,'callback',cb,'backgroundcolor',get(gcf,'color'));

% table ------------------------------------------------------------------------
function h=newtable(tag,pos,n,m,value)
if nargin<5, value=[]; end
w=pos(3)/m; hh=pos(4)/n;
for i=1:n, for j=1:m
    p=[pos(1)+(j-1)*w,pos(2)+(i-1)*hh,w,hh];
    h(i,j)=newedit(sprintf('%s_%d_%d',tag,i,j),p,'',3);
end, end
set(h(1,1),'userdata',[n,m]);
settable(tag,value);

function value=gettable(tag)
nm=getudata([tag,'_1_1']);
for i=1:nm(1), for j=1:nm(2)
    value(i,j)=getnum(sprintf('%s_%d_%d',tag,i,j));
end, end

function settable(tag,value)
nm=getudata([tag,'_1_1']);
for i=1:min(nm(1),size(value,1)), for j=1:min(nm(2),size(value,2))
    setnum(sprintf('%s_%d_%d',tag,i,j),value(i,j));
end, end

% menu -------------------------------------------------------------------------
function h=newmenu(tag,top,labels,vals,cbs)
h=uimenu('tag',tag,'label',top,'handlevisibility','off','visible','off');
addmenus(h,labels,vals,cbs);
updatenow; set(h,'visible','on');

function h=newcmenu(tag,pos,labels,vals,cbs)
h=uicontextmenu('tag',tag,'position',pos,'handlevisibility','off');
addmenus(h,labels,vals,cbs);

function addmenus(h,labels,vals,cbs)
sep='off'; m=1;
for n=1:length(labels)
    if ~strcmp(labels{n},'-')
        if iscell(cbs), cb=cbs{n}; else cb=cbs; end
        if isempty(vals), val=m; else val=vals{n}; end
        if iscell(cb)
            hh=uimenu(h,'label',labels{n},'separator',sep,'userdata',val);
            addmenus(hh,cb{1},cb{2},cb{3}); % submenu
        else
            uimenu(h,'label',labels{n},'callback',cb,'separator',sep,'userdata',val);
        end
        sep='off'; m=m+1;
    else sep='on'; end
end

function showmenu(tag), set(geth(tag),'visible','on');
function hidemenu(tag), set(geth(tag),'visible','off');

function chkmenu(tag,vals)
if nargin<2, allsel=1; else allsel=0; end
for h=get(geth(tag),'children')'
    if allsel|any(strcmp(get(h,'userdata'),vals)), set(h,'checked','on'); end
end

function unchkmenu(tag,vals)
if nargin<2, allsel=1; else allsel=0; end
for h=get(geth(tag),'children')'
    if allsel|any(strcmp(get(h,'userdata'),vals)), set(h,'checked','off'); end
end

function enamenu(tag,vals)
if nargin<2, allsel=1; else allsel=0; end
for h=get(geth(tag),'children')'
    if allsel|any(strcmp(get(h,'userdata'),vals)), set(h,'enable','on'); end
end

function dismenu(tag,vals)
if nargin<2, allsel=1; else allsel=0; end
for h=get(geth(tag),'children')'
    if allsel|any(strcmp(get(h,'userdata'),vals)), set(h,'enable','off'); end
end

function vals=getchk(tag)
vals={};
for h=get(geth(tag),'children')'
    if strcmp(get(h,'checked'),'on'), vals={get(h,'userdata'),vals{:}}; end
end

function on=togglechk(h)
on=strcmp(get(h,'checked'),'off');
if on, set(h,'checked','on'); else set(h,'checked','off'); end

% ok/cancel buttons ------------------------------------------------------------
function h=newokbtn(tag,pos)
h=newbtn(tag,pos,'OK',[mfilename,' cbOk']);

function h=newcancelbtn(tag,pos)
h=newbtn(tag,pos,' Cancel ',[mfilename,' cbCancel']);

function h=newokcancelbtn(tag,pos)
h=newbtnh(tag,pos,{'OK',' Cancel '},{[mfilename,' cbOk'],[mfilename,' cbCancel']});

% list/popup menu --------------------------------------------------------------
function h=newlist(tag,pos,labels,vals,smode,cb)
if nargin<4, vals={}; end
if nargin<5, smode=2; end
if nargin<6, cb=''; end
h=newsel('listbox',tag,pos,labels,vals,smode,cb);

function h=newlistf(tag,pos,labels,vals,smode,cb)
global p_gut
h=newlist(tag,pos,labels,vals,smode,cb);
set(h,'fontname',p_gut.ffont,'fontsize',p_gut.fsize,'backgroundcolor','w');

function h=newpopm(tag,pos,labels,vals,cb)
if nargin<4, vals={}; end
if nargin<5, cb=''; end
h=newsel('popupmenu',tag,pos,labels,vals,1,cb);

function h=newsel(style,tag,pos,labels,vals,smode,cb)
if smode>1, val=[]; else val=1; end
h=newui(style,tag,pos,'max',smode,'callback',cb,'string',labels,...
        'userdata',vals,'value',val,'backgroundcolor','w');

function setlist(tag,labels,vals)
if nargin<3, vals={}; end
h=geth(tag);
if get(h,'max')>1
    set(h,'string',labels,'userdata',vals,'value',[],'listboxtop',1);
else
    if isempty(labels), labels={''}; end
    set(h,'string',labels,'userdata',vals,'value',1,'listboxtop',1);
end

function [labels,vals]=getlist(tag)
h=geth(tag); labels=get(h,'string'); vals=get(h,'userdata');

function setsel(tag,sel)
h=geth(tag);
vals=get(h,'userdata'); if isempty(vals), vals=get(h,'string'); end, i=[];
for n=1:length(vals)
    if iscell(vals), i(n)=any(strcmp(vals{n},sel)); else i(n)=any(vals(n)==sel); end
end
i=find(i);
if ~isempty(i), set(h,'value',i);
elseif get(h,'max')>1, set(h,'value',[]); else set(h,'value',1); end

function [val,label]=getsel(tag)
h=geth(tag); labels=get(h,'string'); vals=get(h,'userdata'); label={};
if isempty(vals), vals=labels; end
if iscell(vals), val={}; else val=[]; end
if isempty(vals), return, end
m=0;
for n=get(h,'value')
    if iscell(vals), val={val{:},vals{n}}; else val=[val,vals(n)]; end
    m=m+1; label{m}=labels{n};
end
if get(h,'max')<=1&iscell(val)
    if isempty(val), val=''; else val=val{1}; end
end

% progress bar -----------------------------------------------------------------
function h=newprog(tag,pos,color)
global p_gut
if nargin<3, color='r'; end
ec=[0.5,0.5,0.5];
h=axes('tag',tag,'xtick',[],'ytick',[],'units','pixels','xlimmode','manual',...
       'ylimmode','manual','box','off','position',pos,'color','w');
rectangle('parent',h,'position',[0,0,0.001,1],'facecolor',color,'edgecolor',ec);
text('parent',h,'position',[0.5,0.5],'horizontal','center',...
     'fontname',p_gut.ufont,'fontsize',p_gut.usize);
line('parent',h,'xdata',[0,0,1],'ydata',[0,1,1],'color',ec)
line('parent',h,'xdata',[0,1,1],'ydata',[0,0,1],'color','w')
updatenow;

function setprog(tag,prog), setprogh(geth(tag),prog)

function setprogh(h,prog)
if isempty(prog), return, end
h1=findhobj(h,'type','text');
h2=findhobj(h,'type','rectangle');
p=max(0,round(prog*100)); x=get(h2,'position');
if p==round(x(3)*100)&~isempty(get(h1,'string')), return, end
set(h1,'string',sprintf('%d%%',p));
set(h2,'position',[0,0,p/100+0.001,1]);
updatenow;

% date input -------------------------------------------------------------------
function h=newymd(tag,pos,ymd,cb)
if nargin<4, cb=''; end
vs={1990:2030,1:12,1:31}; x=[0,47,82]; w=[47,35,35];
pos(1)=pos(1)+min(pos(3)-117,0);
for n=1:3
    strs={''}; for v=vs{n}, strs={strs{:},sprintf('%2d',v)}; end
    h=newpopm([tag,'_',num2str(n)],[pos(1)+x(n),pos(2),w(n),20],strs,{},cb);
end
setymd(tag,ymd);

function ymd=getymd(tag)
for n=1:3, ymd(n)=get(geth([tag,'_',num2str(n)]),'value')-1; end
if all(ymd>0), ymd(1)=ymd(1)+1989; else ymd=[]; end

function setymd(tag,ymd)
if length(ymd)<3, ymd=[0,0,0]; else ymd(1)=ymd(1)-1989; end
for n=1:3, set(geth([tag,'_',num2str(n)]),'value',ymd(n)+1); end

% date/time input --------------------------------------------------------------
function h=newtime(tag,pos,time,opt)
if nargin<4, opt=0; end
if ~opt, x=[0,0,0,118,144,170]; w=25; else x=[0,0,0,118,151]; w=32; end
pos(1)=pos(1)+min(pos(3)-x(end)-w-12,0);
h=newymd(tag,pos+[0,0,60,0],time,'');
for n=4:length(x)
    newedit([tag,'_',num2str(n)],[pos(1)+x(n),pos(2)-2,w,22],0);
end
newbtn([tag,'_btn'],[pos(1)+x(end)+w,pos(2),11,19],'.',[mfilename,' cbShowTime']);
settime(tag,time);

function time=gettime(tag)
time=getymd(tag); if isempty(time), return, end
for n=4:6, time(n)=getnum([tag,'_',num2str(n)]); end

function settime(tag,time)
setymd(tag,time); time(7)=0;
for n=4:6, setstring([tag,'_',num2str(n)],sprintf('%.0f',time(n))); end

function cbShowTime
tag=get(gcbo,'tag');
time=gettime(tag(1:end-4));
[td,ts]=caltomjd(time);
utc_tai=prm_utc_tai(td+ts/86400,1);
leap=19+utc_tai;
tut=mjdtocal(td,ts+leap);
day=td-caltomjd([1980,1,6]);
week=floor(day/7); dow=day-week*7+1;
gpst=(dow-1)*86400+ts;
doy=td-caltomjd([time(1),1,1])+1;
txt1=sprintf(['%04d/%02d/%02d %02d:%02d:%02.0f GPST\n',...
             '%04d/%02d/%02d %02d:%02d:%02.0f UTC \n'],time,tut);
txt2=sprintf('MJD:\nGPS Week:\nGPS Time:\nDay of Year:\nDay of Week:\nTime of Day:\nLeap Seconds:\n');
txt3=sprintf(' %.4f\n %d\n %.0f s\n %03d\n %d\n %.0f s\n%d s\n',...
             td+ts/86400,week,gpst,doy,dow,ts,-leap);
newdlg('','Time',[180,158],1);
newtext('',[5,115,170,35],txt1,2);
newtext('',[5,7,88,110],txt2,3);
newtext('',[93,7,150,110],txt3,1);
newokbtn('',[147,4,30,20]);
if waitok, close, end

% parameters input -------------------------------------------------------------
function h=newprms(tag,pos,prms)
np=size(prms,1); hh=20; if length(pos)<5, w=80; else w=pos(5); end
for n=1:np
    tagn=[tag,'_',num2str(n)];
    v=pos(2)+pos(4)*(np-n); p=[pos(1)+pos(3)-w,v,w,hh];
    t=newtext(['_',tagn],[pos(1),v-7,pos(3),hh+3],prms{n,2});
    switch prms{n,1}
    case 'n', h(n)=newedit(tagn,p,prms{n,3},3);
    case 'e', h(n)=newedit(tagn,p,prms{n,3},1);
    case 'p', h(n)=newpopm(tagn,p+[0,0,1,0],prms{n,3}{1},prms{n,3}{2},'');
    case 'y', h(n)=newymd(tagn,p+[1,0,0,0],prms{n,3},'');
    case 't', h(n)=newtime(tagn,p,prms{n,3});
    case 's', hs=newtable(tagn,p,1,length(prms{n,3}),prms{n,3}); h(n)=hs(1,1);
    case 'f', h(n)=newfedit(tagn,[pos(1),p(2),pos(3)+3,p(4)],prms{n,3});
    case 'h', h(n)=newfedit(tagn,[pos(1),p(2),pos(3)+3,p(4)],prms{n,3},1);
    case 'd', h(n)=newdedit(tagn,[pos(1),p(2),pos(3)+3,p(4)],prms{n,3});
    case 'b', h(n)=newbtn(tagn,p,'...',prms{n,3});
    case 'c', h(n)=newcolorsel(tagn,p,prms{n,3});
    case 'g', h(n)=newfontsel(tagn,p,prms{n,3});
    otherwise, h(n)=-1; end
    if n==1, set(t,'userdata',prms), end
end

function setprms(tag,varargin)
prms=get(geth(['_',tag,'_1']),'userdata');
for n=1:min(size(prms,1),nargin-1)
    tagn=[tag,'_',num2str(n)];
    switch prms{n,1}
    case 'e', setstring(tagn,varargin{n});
    case 'n', setnum(tagn,varargin{n});
    case 'p', setsel(tagn,varargin{n});
    case 'y', setymd(tagn,varargin{n});
    case 't', settime(tagn,varargin{n});
    case 's', settable(tagn,varargin{n});
    case {'f','h','d'}, setstring(tagn,varargin{n});
    case 'c', setcolorsel(tagn,varargin{n});
    case 'g', setfontsel(tagn,varargin{n});
    end
end

function varargout=getprms(tag)
prms=get(geth(['_',tag,'_1']),'userdata');
for n=1:min(size(prms,1),nargout)
    tagn=[tag,'_',num2str(n)];
    switch prms{n,1}
    case 'e', varargout{n}=getstring(tagn);
    case 'n', varargout{n}=getnum(tagn);
    case 'p', varargout{n}=getsel(tagn);
    case 'y', varargout{n}=getymd(tagn);
    case 't', varargout{n}=gettime(tagn);
    case 's', varargout{n}=gettable(tagn);
    case {'f','h','d'}, varargout{n}=char(getstring(tagn));
    case 'c', varargout{n}=getcolorsel(tagn);
    case 'g', varargout{n}=getfontsel(tagn);
    otherwise, varargout{n}=[]; end
end

% file list --------------------------------------------------------------------
function updatefilelist(tag,dirs,patt,order,ud,width,opt)
if nargin<3,patt='*.*'; end
if nargin<4,order=0; end
if nargin<5,ud=0; end
if nargin<6,width=[200,70,100]; end
if nargin<7,opt=0; end
h=gcf; set(h,'pointer','watch')
if opt, f=getudata(tag); else f=dir(fullfile(dirs,patt)); end
f=f(~strncmp({f.name},'.',1)|strcmp({f.name},'..'));
ext={}; time={};
for n=1:length(f)
    [s,ext{n}]=strtok(f(n).name,'.');
    if isfield(f(n),'datenum')
        dv=datevec(f(n).datenum);
        time{n}=sprintf('%02d/%02d/%02d %02d:%02d',mod(dv(1),100),dv(2:5));
    else
        time{n}=datefmt(f(n).date);
    end
end
name={f.name}; for n=find([f.isdir]), name{n}=['<',name{n},'>']; end
switch order
case 1, [s,j]=sort([f.bytes]);
case 2, [s,j]=sort(time);
case 3, [s,j]=sort(ext);
otherwise, [s,j]=sort(name); end
if ud, j=flipud(j(:)); end
f=f(j); time=time(j); name=name(j);
[fn,fs]=getfont('f');
format=sprintf('%%-%.0fs %%%.0fd %%-%.0fs',round(width/(fs*0.79))-1);
str={}; for n=1:length(f), str{n}=sprintf(format,name{n},f(n).bytes,time{n}); end
figure(h)
set(geth(tag),'string',str,'value',[],'listboxtop',1,'userdata',f)
set(h,'pointer','arrow')

function s=datefmt(date)
ms={'Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec'};
mt={'01','02','03','04','05','06','07','08','09','10','11','12'};
s=[date(10:11),'/',mt{strcmp(date(4:6),ms)},'/',date([1:2,12:17])];

function [files,isdir]=getfilelist(tag,dirs)
files={}; isdir=[]; f=getudata(tag); if ~isstruct(f), return, end
for i=get(geth(tag),'value')
    files={files{:},fullfile(dirs,f(i).name)}; isdir=[isdir,f(i).isdir];
end

% directory list ---------------------------------------------------------------
function updatedirlist(tag,dirs)
persistent pdirs, if isempty(pdirs), pdirs=pwd; end
if isempty(dirs), dirs=pdirs; else pdirs=dirs; end
if dirs(end)==':', dirs=[dirs,'\']; end
if strcmp(dirs,'\')
    parent=''; list=getdrvlist; sel={};
elseif findstr(dirs,':\')==length(dirs)-1
    parent=''; list=getdrvlist; sel=find(strcmp(list,lower(dirs)));
else
    [parent,d,e]=fileparts(dirs);
    f=dir(parent); 
    list={'..',f([f.isdir]&~strncmp({f.name},'.',1)).name};
    sel=find(strcmp(lower(list),lower([d,e])));
end
if isempty(sel), sel=1; end
set(geth(tag),'string',list,'value',sel,'userdata',parent)

function dirs=getdirlist(tag)
parent=get(geth(tag),'userdata');
list=get(geth(tag),'string');
if isempty(list), dirs=pwd;
else dirs=fullfile(parent,list{get(geth(tag),'value')}); end

function ds=getdrvlist
persistent drv, if ~isempty(drv), ds=drv; return; end
drv={'a:\'};
for d=double('b':'z')
    if ~isempty(dir([char(d),':\'])), drv={drv{:},[char(d),':\']}; end
end
ds=drv;

% directory input edit ---------------------------------------------------------
function h=newdedit(tag,pos,str)
h=newedit(tag,pos-[0,0,19,0],str,1);
b=newbtn([tag,'b'],[pos(1)+pos(3)-19,pos(2)+1,15,pos(4)-2],'...',[mfilename,' cbDirEdit']);
set(b,'userdata',h,'fontsize',8)

function cbDirEdit
[ok,dirs]=dirdlg(get(get(gcbo,'userdata'),'string'));
if ok, set(get(gcbo,'userdata'),'string',dirs), end

% file input edit -------------------------------------------------------------
function h=newfedit(tag,pos,str,opt)
if nargin<4, opt=0; end
if opt
    h=newedit(tag,pos-[0,0,36,0],str,1);
    b=newbtnh([tag,'b'],[pos(1)+pos(3)-37,pos(2)+1,34,pos(4)-2],{'...','.'},...
              {[mfilename,' cbFileSel'],[mfilename,' cbFileEdit']});
else
    h=newedit(tag,pos-[0,0,19,0],str,1);
    b=newbtn([tag,'b'],[pos(1)+pos(3)-19,pos(2)+1,15,pos(4)-2],'...',...
             [mfilename,' cbFileSel']);
end
set(b,'userdata',h,'fontsize',8)

function cbFileSel
file=gfilepath('',get(get(gcbo,'userdata'),'string'));
if ~isempty(file), [d,f]=fileparts(file); wd=cd; if ~isempty(d), cd(d); end, end
[f,path]=uigetfile('*.*','Select File');
if f~=0, set(get(gcbo,'userdata'),'string',fullfile(path,f)), end
if ~isempty(file), cd(wd); end

function cbFileEdit
global p_gut
file=gfilepath('',get(get(gcbo,'userdata'),'string'));
if ~exist(file), return, end
[s,w]=dos(['"',p_gut.editor,'" "',file,'" &; exit']);
if s, errdlg(w); end

% font select edit ------------------------------------------------------------
function h=newfontsel(tag,pos,font)
h1=newtext('',[0,0,1,1],'');
h2=newtext('',[pos(1)+pos(3)-182,pos(2)-3,160,pos(4)],'',3);
h=newbtn(tag,[pos(1)+pos(3)-16,pos(2)+2,15,18],'...',[mfilename,' cbSelFont']);
set(h,'userdata',[h1,h2],'fontsize',8)
if isstruct(font), setfontselh(h,font); end

function cbSelFont, setfontselh(gcbo,uisetfont(getfontselh(gcbo)));

function setfontselh(h,font)
if ~isstruct(font), return, end
h=get(h,'userdata');
set(h(1),font);
set(h(2),'string',sprintf('%s %dpt',font.FontName,font.FontSize),...
    'fontname',font.FontName,'fontsize',font.FontSize,...
    'fontangle',font.FontAngle,'fontweight',font.FontWeight);

function font=getfontselh(h), h=get(h,'userdata'); font=getfonth(h(1));

function setfontsel(tag,font), setfontselh(geth(tag),font)
function font=getfontsel(tag), font=getfontselh(geth(tag));

% color select edit -----------------------------------------------------------
function h=newcolorsel(tag,pos,color)
h=newtoggle([tag,'_2'],[pos(1),pos(2),14,pos(4)],[mfilename,' cbEnaColor']);
if isnan(color), off=0; set(h,'visible','off'); else off=15; end
h=newbtn([tag,'_1'],[pos(1)+pos(3)-16,pos(2),16,pos(4)],'...',[mfilename,' cbSelColor']);
set(h,'fontsize',8);
h=newfrm([tag,'_3'],pos+[off,1,-off-17,-2]);
setcolorsel(tag,color);

function cbSelColor
tag=get(gcbo,'tag');
setcolorsel(tag(1:end-2),uisetcolor(getcolorsel(tag(1:end-2)),'Select Color'));

function cbEnaColor
tag=get(gcbo,'tag');
if ~get(gcbo,'value'), color='none';
else c=get(gcbo,'userdata'); if isempty(c), color=[0.7,0.7,0.7]; else color=c; end, end
setcolorsel(tag(1:end-2),color)

function setcolorsel(tag,color)
if isempty(color)|isnan(color)|strcmp(color,'none')
    set(geth([tag,'_1']),'enable','off');
    set(geth([tag,'_3']),'foregroundcolor',[0.8,0.8,0.8],'backgroundcolor',get(gcf,'color'));
    setval([tag,'_2'],0);
else
    set(geth([tag,'_1']),'enable','on');
    set(geth([tag,'_3']),'foregroundcolor','k','backgroundcolor',color);
    set(geth([tag,'_2']),'value',1,'userdata',color);
end

function color=getcolorsel(tag)
if ~get(geth([tag,'_2']),'value'), color='none';
else color=get(geth([tag,'_3']),'backgroundcolor'); end

% text viewer -----------------------------------------------------------------
function h=newviewer(tag,title,siz,varargin)
global p_gut
str={}; file={}; callback=''; n=1;
while n<=nargin-3
    switch varargin{n}
    case 'str',  str =varargin{n+1}; n=n+2;
    case 'file', file=varargin{n+1}; n=n+2;
    case 'callback', callback=varargin{n+1}; n=n+2;
    otherwise, n=n+1; end
end
h=newfig(tag,title,siz,[mfilename,' cbViewClose'],'position',p_gut.viewpos);
if ~isempty(h)
    set(h,'resize','on','ResizeFcn',[mfilename,' cbViewResize']);
    file=unique(file);
    if length(file)>1
        newpopm('vfile',[0,0,100,23],file,{},[mfilename,' cbViewFile']);
        setsel('vfile',file{end});
        newbtnh('btnv',[0,0,70,20],{'<','>'},{[mfilename,' cbViewFB'],[mfilename,' cbViewFF']});
    end
    newlistf('vstr',[0,0,300,300],{},{},2,callback);
    newbtnh('btns',[0,0,210,20],{' Edit... ',' Save... ',' Option... '},...
       {[mfilename,' cbViewEdit'],[mfilename,' cbViewSave'],[mfilename,' cbViewOpt']});
    newbtn('btnf',[0,0,40,20],' Find ',[mfilename,' cbViewFind']);
    newbtn('btnc',[0,0,70,20],'Close','');
    newedit('fstr',[0,0,100,21],'',1);
    if length(file)==0, setdis('btns_1'); end
else
    set(geth(tag),'name',title);
    if length(file)>1, setlist('vfile',file,{}); end
end
resizeviewer(h);
if p_gut.viewln, for n=1:length(str), str{n}=[sprintf('%5d: ',n),str{n}]; end, end
set(geth('vstr'),'string',str,'backgroundcolor',p_gut.viewbc,...
    'foregroundcolor',p_gut.viewfc,p_gut.viewfont);
if length(file)>1, cbViewFile, elseif length(file)==1, ViewFile(file{1}); end

function resizeviewer(h)
p=get(h,'position'); p(3:4)=max(p(3:4),[100,100]);
for n=1:3, setpos(['btns_',num2str(n)],[(n-1)*70+2,p(4)-22]); end
setpos('fstr',[215,p(4)-23]);
setpos('btnf',[316,p(4)-22]);
setpos('btnc',[p(3)-71,p(4)-22]); p(4)=p(4)-24;
h=geth('vfile');
if ~isempty(h)
    for n=1:2, setpos(['btnv_',num2str(n)],[(n-1)*35+365,p(4)+2]); end
    set(h,'position',[1,p(4)-23,p(3),23]);
    p(4)=p(4)-22;
end
set(geth('vstr'),'position',[0,0,p(3)+2,p(4)]);

function cbViewResize
resizeviewer(gcbo);

function cbViewClose
global p_gut
p_gut.viewpos=get(gcf,'position');
closereq

function cbViewFB
h=geth('vfile');
set(h,'value',max(get(h,'value')-1,1)); cbViewFile

function cbViewFF
h=geth('vfile'); n=length(get(h,'string'));
set(h,'value',min(get(h,'value')+1,n)); cbViewFile

function cbViewOpt
global p_gut
h=geth('vstr');
prm1={
'c','Background Color', nan
'c','Foreground Color', nan
'g','Text Font',        []
'p','Show Line Number', {{'OFF','ON'},0:1}
};
newdlg('','Text Viewer : Options',[250,147],1);
newprms('prm1',[13,52,232,23],prm1);
newtext('',[5,24,240,20],['(Matlab : ',version,')'],2);
setprms('prm1',p_gut.viewbc,p_gut.viewfc,p_gut.viewfont,p_gut.viewln);
gut('newokcancelbtn','',[65,4,180,22]);
if ~waitok, return, end
[p_gut.viewbc,p_gut.viewfc,p_gut.viewfont,viewln]=getprms('prm1');
close
set(h,'backgroundcolor',p_gut.viewbc,'foregroundcolor',p_gut.viewfc,p_gut.viewfont);
if viewln==p_gut.viewln, return, end
str=get(h,'string');
for n=1:length(str)
    if viewln, str{n}=[sprintf('%5d: ',n),str{n}]; else str{n}=str{n}(8:end); end
end
set(h,'string',str);
p_gut.viewln=viewln;

function cbViewFile
ViewFile(getsel('vfile'));
h=geth('vfile'); p=get(h,'value'); n=length(get(h,'string'));
if p==1, setdis('btnv_1'); else setena('btnv_1'); end
if p==n, setdis('btnv_2'); else setena('btnv_2'); end

function cbViewSave
global p_gut
persistent ppath
if ~isempty(ppath), wd=cd; cd(ppath); else wd=''; end
[file,path]=uiputfile('*.*','Save To Text File');
if ~isempty(wd), cd(wd); end
if file==0, return, else ppath=path; end
f=fopen([path,file],'wt'); if f==-1, return, end
for s=getstring('vstr')'
    if p_gut.viewln, str=s{1}(8:end); else str=s{1}; end
    fprintf(f,'%s\n',str);
end
fclose(f);

function cbViewEdit
global p_gut
[t,file]=strtok(get(gcf,'name'),':');
[s,w]=dos(['"',p_gut.editor,'" "',file(3:end),'" &; exit']);
if s, errdlg(w); end

function ViewFile(file)
global p_gut
file=gfilepath('',file);
if ~exist(file), return, end
h=gcf; set(h,'pointer','watch')
try
    org=file; [file,stat,msg]=uncompact(org,'org','dir',pwd);
    if stat==0
        str=textread(file,'%s','delimiter','\n','whitespace','','bufsize',256);
        if length(str)>p_gut.maxtext
            str={str{1:p_gut.maxtext},'...',...
                 ['(truncated in ',num2str(p_gut.maxtext),' lines)']};
        end
        if p_gut.viewln
            for n=1:length(str), str{n}=[sprintf('%5d: ',n),str{n}]; end
        end
        set(geth('vstr'),'string',str,'value',[]);
        set(gcf,'name',[strtok(get(gcf,'name'),':'),': ',org])
        
        if ~strcmp(file,org), gut('setdis','btns_1'); delete(file); end
    else
        errdlg([msg,' : ',file]);
    end
catch
    err=lasterror; errdlg([err.message,' : ',file]);
end
set(h,'pointer','arrow')

function cbViewFind
fstr=getstring('fstr');
vstr=getlist('vstr');
pos=get(geth('vstr'),'value');
if isempty(pos)|length(vstr)<=pos, pos=1; else pos=pos+1; end
for n=pos:length(vstr)
    if ~isempty(findstr(fstr,vstr{n})), set(geth('vstr'),'value',n); return, end
end
set(geth('vstr'),'value',[]);

% message and progress bar -----------------------------------------------------
function h=newmsgbar(tag,title,siz,color,abtn)
if nargin<4, color=[0,0.7,1]; end
if nargin<5, abtn=0; end
h=newdlg(tag,title,siz);
set(h,'doublebuffer','on');
newtext('msg',[5,siz(2)/2-2,siz(1)-10,18],'',2);
newprog('prog',[10,7,siz(1)-20,15],color);
if abtn
    hh=newbtn('abtn',[siz(1)-63,4,60,20],'Abort',[mfilename,' cbMsgAbort']);
    set(hh,'userdata',0);
    setsize('prog',[siz(1)-78,15])
end

function abort=setmsgbar(tag,msg,prog)
h=findhobj(0,'tag',tag); if isempty(h), abort=1; return, end
hh=findhobj(h,'tag','msg'); if ~isnan(msg), set(hh,'string',msg); end
hh=findhobj(h,'tag','abtn'); abort=get(hh,'userdata');
hh=findhobj(h,'tag','prog'); setprogh(hh,prog);

function closemsgbar(tag)
h=findhobj(0,'tag',tag); if isempty(h), return, end
delete(h);

function cbMsgAbort, set(gcbo,'userdata',1,'enable','off');

% select directory dialog ------------------------------------------------------
function [ok,dirs]=dirdlg(dirs)
persistent hist, if isempty(hist), hist={''}; end
if isempty(hist{1}), dirs=pwd; else dirs=hist{1}; end
newdlg('','Select Directory',[285,300]);
newbtnh('',[4,277,48,18],{'\','<','>'},{[mfilename,' cbDirRoot'],...
        [mfilename,' cbDirUp'],[mfilename,' cbDirDown']});
newcmenu('drvs',[4,277],getdrvlist,{},[mfilename,' cbDirDrv']);
newpopm('path',[50,277,235,20],{dirs,hist{:}},{},[mfilename,' cbDirSel']);
newlist('dirs',[2,28,284,246],{},{},1,'cbDirDlgClk');
newbtn('',[5,4,92,22],'New...',[mfilename,' cbMkDir']);
newokcancelbtn('',[98,4,184,22]);
UpdateDirDlg(dirs)
ok=waitok; if ~ok, dirs=''; return, end
dirs=getdirlist('dirs');
[p,f,e]=fileparts(dirs); if strcmp([f,e],'..'), dirs=p; end
if ~any(strcmp(dirs,hist)), hist={dirs,hist{:}}; end
close

function cbDirDlgClk
dirs=getdirlist('dirs');
if strcmp(get(gcf,'selectiontype'),'open')
    [parent,f,e]=fileparts(dirs);
    if strcmp([f,e],'..'), dirs=parent; else dirs=fullfile(dirs,'..'); end
end
UpdateDirDlg(dirs)

function cbDirRoot
setvis('drvs');

function cbDirDrv
UpdateDirDlg([get(gcbo,'label'),'..']);

function cbDirUp
dirs=gut('getdirlist','dirs');
[parent,file,ext]=fileparts(dirs);
UpdateDirDlg(parent);

function cbDirDown
dirs=gut('getdirlist','dirs');
[p,f,e]=fileparts(dirs);
if ~strcmp([f,e],'..'), UpdateDirDlg(fullfile(dirs,'..')); end

function cbDirSel
dirs=getsel('path'); if exist(dirs), UpdateDirDlg(dirs); end

function cbMkDir
ndir=newdirdlg; if isempty(ndir), return, end
[parent,f,e]=fileparts(getdirlist('dirs'));
cdir=pwd; cd(parent),
if mkdir(ndir)==1, UpdateDirDlg(fullfile(parent,ndir)), end
cd(cdir)

function UpdateDirDlg(dirs)
updatedirlist('dirs',dirs)
dirs=getdirlist('dirs');
list=getlist('path'); list{1}=dirs;
setlist('path',list)

% new directory dialog ---------------------------------------------------------
function ndir=newdirdlg
newdlg('','New Directory',[240,55]);
h=newedit('ndir',[5,30,230,20],'',1);
set(h,'callback',[mfilename,' cbOk'])
newokcancelbtn('',[75,5,160,20]);
if ~waitok, ndir=''; return, end
ndir=getstring('ndir');
close

% about dialog -----------------------------------------------------------------
function aboutdlg(title,siz,msgs)
newdlg('',title,siz,0);
h=newtext('',[0,siz(2)-45,siz(1),23],msgs{1},2);
set(h,'fontsize',10);
for n=2:length(msgs)
    h=newtext('',[0,siz(2)-32-n*18,siz(1),20],msgs{n},2);
end
newokbtn('',[(siz(1)-80)/2,8,80,20]);
if waitok, close, end

% confirmation dialog ----------------------------------------------------------
function ok=confdlg(msg,btns)
if nargin<2, btns={'OK',' Cancel '}; end
newdlg('','',[320,80]);
newtext('',[10,35,320,40],msg,2);
newbtnh('',[70,4,180,22],btns,{[mfilename,' cbOk'],[mfilename,' cbCancel']});
ok=waitok; if ok, close, end

% error dialog -----------------------------------------------------------------
function ok=errdlg(msg)
newdlg('','Error',[340,92]);
newtext('',[10,30,320,56],msg,2);
newokbtn('',[130,4,80,22]);
ok=waitok; if ok, close, end

% select modal dialog ----------------------------------------------------------
function [sel,ok]=seldlg(title,siz,labels,vals,sel,smode)
if nargin<4, vals={}; end
if nargin<5, sel={};  end
if nargin<6, smode=1; end
newdlg('',title,siz);
newlist('list',[2,28,siz(1)-2,siz(2)-29],labels,vals,smode,'cbSelDlgClk');
setsel('list',sel);
newokcancelbtn('',[siz(1)-185,3,180,22]);
ok=waitok; if ~ok, sel={}; return, end
sel=getsel('list');
close

function cbSelDlgClk, if strcmp(get(gcf,'selectiontype'),'open'), cbOk, end

% range input dialog -----------------------------------------------------------
function [range,ok]=rangedlg(title,siz,range)
newdlg('',title,siz);
y=(siz(2)-25)/2+15; w=(siz(1)-36)/2;
newtext('',[siz(1)/2-10,y,20,20],'-',2);
newedit('r1',[10,  y,w,22],'');
newedit('r2',[27+w,y,w,22],'');
newokcancelbtn('',[siz(1)-165,3,160,21]);
if ~isempty(range), setnum('r1',range(1)); setnum('r2',range(2)); end
ok=waitok; if ~ok, return, end
range=[getnum('r1'),getnum('r2')];
close

% wait for push button ---------------------------------------------------------
function ok=waitok
h=gcf; waitfor(h,'userdata'), if ~any(get(0,'children')==h), ok=0; return, end
ok=get(h,'userdata'); if ~ok, close(h), end

% common callback --------------------------------------------------------------
function cbOk, set(gcf,'userdata',1)
function cbCancel, set(gcf,'userdata',0)

% update now -------------------------------------------------------------------
function updatenow
h=gcf; drawnow; figure(h);

% string<->number --------------------------------------------------------------
function num=StrToNum(str)
if isempty(str), num=0; else num=str2num(str); if isempty(num), num=0; end, end

function str=NumToStr(num), str=num2str(num,6);
