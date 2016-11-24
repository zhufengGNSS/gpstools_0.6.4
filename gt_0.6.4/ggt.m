function varargout=ggt(varargin)
%-------------------------------------------------------------------------------
% [system] : GpsTools
% [module] : graph common libraries for GpsTools
% [func]   : graph common libraries for GpsTools
% [argin]  : function,argins
%              'newaxes',tag,pos(,opts)     : generate axes
%              'setaxes',tag,name,val(,...) : set axes properties
%              'getaxes',tag,name           : get axes properties
%              'subplotv',tag,n,np(,opts)   : generate suplot axpes vertical
%              'subploth',tag,n,np(,opts)   : generate suplot axes horizontal
%              'subplot',tag,n,m,np,mp(,opts) : generate suplot axes
%              'stext',tag,pos,str,align(,opts) : static text
%              'mtext',tag,str,align(,opts) : add moving text
%              'delmtext',tag               : delete moving text
%              'skymap'(,opts)              : draw skymap
%              'skyplot',azel(,opts)        : plot in skymap
%              'skypatch',azel,c,(,opts)    : patch in skymap
%              'skytext',azel,str,(,opts)   : text in skymap
%              'skycpatch',azel,wa,we,c,(,opts) : corn patch in skymap
%              'colorbarv',pos,range,label(,opts) : colorbar vertical
%              'colorbarh',pos,range,label(,opts) : colorbar horizontal
%              'barplot',x,labels(,opts)    : plot in bargraph
%              'cmapdlg',tag,title          : colormap dialog
% [argout] : argouts
% [note]   :
% [version]: $Revision: 12 $ $Date: 2008-11-25 10:02:15 +0900 (火, 25 11 2008) $
%            Copyright(c) 2004-2006 by T.Takasu, all rights reserved
% [history]: 04/12/17  0.1  new
%            05/06/28  0.2  expand taxis range
%            06/03/06  0.3  add skymap opts dirs, rotate
%-------------------------------------------------------------------------------
if nargout==0, feval(varargin{:});
else [varargout{1:nargout}]=feval(varargin{:}); end

% axes -------------------------------------------------------------------------
function ax=newaxes(tag,pos,varargin)
[fname,fsize]=gut('getfont','g');
prm.equal=0;             % axis equal flag
prm.move =[0,0,0,0];     % moving axis flag [x,y,xzoom,yzoom]
prm.moved=0;             % moving axis direction (0:no move)
prm.movel=[0,0,0,0];     % moving axis link [x,y,xzoom,yzoom]
prm.point=[0,0];         % base point
prm.limit=[0,0,0,0];     % base axis limit
prm.td   =[];            % time axis base date (mjd)
prm.topts={};            % time axis option
prm.texth=[];            % moving text handles
prm.texta=[];            % moving text alignments
prm.callback='';         % callback on move end
prm.callbacks='';        % callback on move start
prm.callbackm='';        % callback on move
prm.callbackd='';        % callback on double click
va={}; n=1;
while n<=nargin-2
    if ischar(varargin{n})
        switch varargin{n}
        case 'equal', prm.equal=1; n=n+1;
        case 'move',  prm.move =varargin{n+1}; n=n+2;
        case 'link',  prm.movel=varargin{n+1}; n=n+2;
        case 'taxis', prm.td   =varargin{n+1}; n=n+2;
        case 'topts', prm.topts=varargin{n+1}; n=n+2;
        case 'callback', prm.callback=varargin{n+1}; n=n+2;
        case 'callbacks', prm.callbacks=varargin{n+1}; n=n+2;
        case 'callbackm', prm.callbackm=varargin{n+1}; n=n+2;
        case 'callbackd', prm.callbackd=varargin{n+1}; n=n+2;
        case {'xlim','ylim'},
            if varargin{n+1}(1)>=varargin{n+1}(2)
                varargin{n+1}(2)=varargin{n+1}(1)+1;
            end
            va={va{:},varargin{n:n+1}}; n=n+2;
        otherwise va={va{:},varargin{n}}; n=n+1; end
    else
        va={va{:},varargin{n}}; n=n+1;
    end
end
ax=axes('position',pos); hold on, grid on, box on
if prm.equal, axis equal, end
set(ax,'tag',['ggt_',tag],'fontname',fname,'fontsize',fsize,'userdata',prm,va{:});
UpdateTaxis(ax);

if any(prm.move)
    set(gcf,'windowbuttondownfcn',  [mfilename,' cbBtnDown'],...
            'windowbuttonmotionfcn',[mfilename,' cbBtnMove'],...
            'windowbuttonupfcn',    [mfilename,' cbBtnUp'  ]);
end

% get/set axes properties ------------------------------------------------------
function setaxes(tag,varargin)
set(geth(['ggt_',tag]),varargin{:});

function value=getaxes(tag,name)
value=get(geth(['ggt_',tag]),name);

% subplot vertical/horizontal --------------------------------------------------
function ax=subplot(tag,nn,mm,np,mp,varargin)
margin=[0.1,0.01,0.1,0.01]; % [fig l/r,axis l/r,fig t/b,axis,t/b]
va={}; n=1;
while n<=nargin-5
    if ischar(varargin{n})
        switch varargin{n}
        case 'margin', margin=varargin{n+1}; n=n+2;
        otherwise va={va{:},varargin{n}}; n=n+1; end
    else
        va={va{:},varargin{n}}; n=n+1;
    end
end
w=(1-2*margin(1))/np;
h=(1-2*margin(3))/mp;
pos=[(nn-1)*w+margin(1)+margin(2),(mp-mm)*h+margin(3)+margin(4),...
     w-2*margin(2),h-margin(4)];
ax=newaxes([tag,'_',num2str(nn),'_',num2str(mm)],pos,va{:});

function ax=subplotv(tag,nn,np,varargin)
margin=[0.1,0.07,0.05,0.03]; % [left,right,figure top/bottom,axes top/bottom]
va={}; n=1;
while n<=nargin-3
    if ischar(varargin{n})
        switch varargin{n}
        case 'margin', margin=varargin{n+1}; n=n+2;
        otherwise va={va{:},varargin{n}}; n=n+1; end
    else
        va={va{:},varargin{n}}; n=n+1;
    end
end
h=(1-2*margin(3))/np;
pos=[margin(1),(np-nn)*h+margin(3)+margin(4),1-margin(1)-margin(2),h-2*margin(4)];
ax=newaxes([tag,'_',num2str(nn)],pos,va{:});

function ax=subploth(tag,nn,np,varargin)
margin=[0.1,0.07,0.05,0.03]; % [bottom,top,figure left/right,axes left/right]
va={}; n=1;
while n<=nargin-3
    if ischar(varargin{n})
        switch varargin{n}
        case 'margin', margin=varargin{n+1}; n=n+2;
        otherwise va={va{:},varargin{n}}; n=n+1; end
    else
        va={va{:},varargin{n}}; n=n+1;
    end
end
w=(1-2*margin(3))/np;
pos=[(nn-1)*w+margin(3)+margin(4),margin(1),w-2*margin(4),1-margin(1)-margin(2)];
ax=newaxes([tag,'_',num2str(nn)],pos,va{:});

% static text ------------------------------------------------------------------
function h=stext(tag,pos,str,align,varargin)
if nargin<4, align=1; end  % algn=1:left,2:right,3:top center,4:bottom center
ax=gut('geth',['ggt_',tag]); if isempty(ax), return, end
fname=get(ax,'fontname');
fsize=get(ax,'fontsize');
punit=GetPixUnit(ax);
marginx=punit(1)*0;
marginy=punit(2)*0;
switch align
case 1, mx= marginx; my=0; va='middle'; ha='left';
case 2, mx=-marginx; my=0; va='middle'; ha='right';
case 3, mx=0; my=-marginy; va='top';    ha='center';
case 4, mx=0; my= marginy; va='bottom'; ha='center';
end
h=text(pos(1)+mx,pos(2)+my,str,'horizontal',ha,'vertical',va,'fontname',fname,...
       'fontsize',fsize,varargin{:});

% moving text ------------------------------------------------------------------
function h=mtext(tag,str,align,varargin)
if nargin<3, align=1; end
ax=gut('geth',['ggt_',tag]); if isempty(ax), return, end
fname=get(ax,'fontname');
fsize=get(ax,'fontsize');
[pos,ha,va]=TextPos(ax,align);
prm=get(ax,'userdata');
h=text(pos(1),pos(2),str,'horizontal',ha,'vertical',va,'fontname',fname,...
       'fontsize',fsize,'interpreter','none',varargin{:});
prm.texth=[prm.texth,h];
prm.texta=[prm.texta,align];
set(ax,'userdata',prm);

% delete mtext -----------------------------------------------------------------
function delmtext(tag)
ax=gut('geth',['ggt_',tag]); if isempty(ax), return, end
prm=get(ax,'userdata'); if isempty(prm.texth), return, end
for h=prm.texth, delete(h); end
prm.texth=[];
prm.texta=[];
set(ax,'userdata',prm);

% callback on button down ------------------------------------------------------
function cbBtnDown
for ax=get(gcf,'children')'
    if ishandle(ax)&findstr(get(ax,'tag'),'ggt_')
        if strcmp(get(gcf,'selectiontype'),'open') % double click
            prm=get(ax,'userdata');
            if ~isempty(prm.callbackd), eval(prm.callbackd); end
        else
            StartMoveAxes(ax);
        end
    end
end

% callback on button move ------------------------------------------------------
function cbBtnMove
for ax=get(gcf,'children')'
    if ishandle(ax)&findstr(get(ax,'tag'),'ggt_'), MoveAxes(ax); end
end

% callback on button up --------------------------------------------------------
function cbBtnUp
for ax=get(gcf,'children')'
    if ishandle(ax)&findstr(get(ax,'tag'),'ggt_'), StopMoveAxes(ax); end
end

% start move axes --------------------------------------------------------------
function StartMoveAxes(ax)
prm=get(ax,'userdata'); if all(~prm.move), return, end
on=IsOnAxis(ax); if on==0|on>=3, return, end
if strcmp(get(gcf,'selectiontype'),'alt'), moved=on+2; else moved=on; end
if ~prm.move(moved), return, end
prm.moved=moved;
prm.point=GetPoint(ax);
prm.limit=[get(ax,'xlim'),get(ax,'ylim')];
set(ax,'userdata',prm);
eval(prm.callbacks);
SetPointer(moved)

% move axes --------------------------------------------------------------------
function MoveAxes(ax)
prm=get(ax,'userdata'); if prm.moved<=0, return, end
dpoint=GetPoint(ax)-prm.point;
fp=0.1.^dpoint; fp=min(max(fp,0.1),10);
xl=prm.limit(1:2);
yl=prm.limit(3:4);
switch prm.moved
case 1, xl=xl-(xl(2)-xl(1))*dpoint(1);
case 2, yl=yl-(yl(2)-yl(1))*dpoint(2);
case 3, xl=Zoom(xl,0.5,fp(1)); if prm.equal, yl=Zoom(yl,0.5,fp(1)); end
case 4, yl=Zoom(yl,0.5,fp(2)); if prm.equal, xl=Zoom(xl,0.5,fp(2)); end
end
UpdateAx(ax,xl);
UpdateAy(ax,yl);
if prm.movel(prm.moved)
    for axx=get(gcf,'children')' % link to other axes
        if ishandle(axx)&findstr(get(axx,'tag'),'ggt_')
            prmx=get(axx,'userdata');
            if prmx.move(prm.moved)
                if mod(prm.moved,2), UpdateAx(axx,xl); else UpdateAy(axx,yl); end
            end
        end
    end
end
eval(prm.callbackm);

function UpdateAx(ax,xl)
if all(~isinf(xl)), set(ax,'xlim',xl); end
UpdateTaxis(ax);
UpdateMtext(ax);

function UpdateAy(ax,yl)
if all(~isinf(yl)), set(ax,'ylim',yl); end
UpdateMtext(ax);

function lim=Zoom(lim,cent,fact)
lim=lim(1)+(lim(2)-lim(1))*(cent-[cent,cent-1]*fact);

% stop move axes ---------------------------------------------------------------
function StopMoveAxes(ax)
prm=get(ax,'userdata'); if prm.moved<=0, return, end
set(gcf,'pointer','arrow')
eval(prm.callback);
prm.moved=0;
set(ax,'userdata',prm);
set(gcf,'pointer','arrow')

% test point on axis -----------------------------------------------------------
function on=IsOnAxis(ax)
margin=GetPixUnit(ax)*10;
xl=get(ax,'xlim');
yl=get(ax,'ylim');
p=get(ax,'currentpoint'); p=p(1,1:2);
on=[1,2,3,4];
if strcmp(get(ax,'xaxislocation'),'top'  ), on([1,3])=[3,1]; end
if strcmp(get(ax,'yaxislocation'),'right'), on([2,4])=[4,2]; end
if     abs(p(2)-yl(1))<margin(2)&xl(1)<=p(1)&p(1)<=xl(2), on=on(1); % bottom
elseif abs(p(1)-xl(1))<margin(1)&yl(1)<=p(2)&p(2)<=yl(2), on=on(2); % left
elseif abs(p(2)-yl(2))<margin(2)&xl(1)<=p(1)&p(1)<=xl(2), on=on(3); % top
elseif abs(p(1)-xl(2))<margin(1)&yl(1)<=p(2)&p(2)<=yl(2), on=on(4); % right
else on=0; end

% mtext position in axes -------------------------------------------------------
function [pos,ha,va]=TextPos(ax,align)
punit=GetPixUnit(ax);
marginx=5*punit(1);
marginy=3*punit(2);
xl=get(ax,'xlim');
yl=get(ax,'ylim');
pos=[]; ha=''; va='';
switch align
case 1, pos=[xl(2)-marginx,yl(2)-marginy]; ha='right'; va='top';    % right top
case 2, pos=[xl(1)+marginx,yl(2)-marginy]; ha='left';  va='top';    % left top
case 3, pos=[xl(2)-marginx,yl(1)+marginy]; ha='right'; va='bottom'; % right bottom
case 4, pos=[xl(1)+marginx,yl(1)+marginy]; ha='left';  va='bottom'; % left bottom
case 5, pos=[xl(2)-marginx,yl(2)]; ha='right'; va='bottom';         % right upper
case 6, pos=[xl(1)+marginx,yl(2)]; ha='left';  va='bottom';         % left upper
end

% get current point in normalized axes coordinate ------------------------------
function p=GetPoint(ax)
pos=get(gcf,'position');
p=get(gcf,'currentpoint')./pos(3:4);
pos=get(ax,'position');
p=(p-pos(1:2))./pos(3:4);

% get axes pixel units ---------------------------------------------------------
function punit=GetPixUnit(ax)
units=get(ax,'units');
set(ax,'units','pixel');
siz=get(ax,'position');
set(ax,'units',units);
xl=get(ax,'xlim');
yl=get(ax,'ylim');
punit=[(xl(2)-xl(1))/siz(3),(yl(2)-yl(1))/siz(4)];

% set pointer ------------------------------------------------------------------
function SetPointer(moved)
if     moved==1, set(gcf,'pointer','right'), return,
elseif moved==2, set(gcf,'pointer','top'),   return, end
o=nan;
p=[
o,o,o,o,o,o,o,2,1,2,o,o,o,o,o,o
o,o,o,o,o,o,2,1,1,1,2,o,o,o,o,o
o,o,o,o,o,2,1,1,1,1,1,2,o,o,o,o
o,o,o,o,2,1,1,1,1,1,1,1,2,o,o,o
o,o,o,2,1,1,1,1,1,1,1,1,1,2,o,o
o,o,o,2,2,2,2,2,1,2,2,2,2,2,o,o
o,o,o,o,o,o,o,2,1,2,o,o,o,o,o,o
o,o,o,o,o,o,o,2,1,2,o,o,o,o,o,o
o,o,o,o,o,2,2,2,1,2,2,2,o,o,o,o
o,o,o,o,o,2,1,1,1,1,1,2,o,o,o,o
o,o,o,o,o,2,2,2,1,2,2,2,o,o,o,o
o,o,o,o,o,o,o,2,1,2,o,o,o,o,o,o
o,o,o,o,o,2,2,2,1,2,2,2,o,o,o,o
o,o,o,o,o,2,1,1,1,1,1,2,o,o,o,o
o,o,o,o,o,o,2,1,1,1,2,o,o,o,o,o
o,o,o,o,o,o,o,2,1,2,o,o,o,o,o,o
];
if moved==3, p=fliplr(p'); end
set(gcf,'pointer','custom','pointershapecdata',p,'pointershapehotspot',[9,9])

% update moving text -----------------------------------------------------------
function UpdateMtext(ax)
prm=get(ax,'userdata');
for n=1:length(prm.texth)
    set(prm.texth(n),'position',TextPos(ax,prm.texta(n)));
end

% update time axes -------------------------------------------------------------
function UpdateTaxis(ax)
prm=get(ax,'userdata'); if isempty(prm.td), return, end
ttt=[3240,1080,360,120,96,48,24,12,6,3,1,1/2,1/6,1/12,1/30,1/60,1/120,1/360,1/720,1/1800,1/3600];
tspan=get(ax,'xlim');
i=max(find(ttt>=(tspan(2)-tspan(1))/32));
if isempty(i), tt=6480; else tt=ttt(i); end
xt=(floor(tspan(1)/tt):floor(tspan(2)/tt))*tt;
xl={};
if all(~strcmp(prm.topts,'nolabel'))
    nx=length(xt); tk=tt; 
    if tt>=120
        if nx>32, tk=tt*8; elseif nx>16, tk=tt*4; elseif nx>8, tk=tt*2; end
    elseif tt>=12
        if nx>16, tk=tk*2; end
    else
        if nx>16, tk=ttt(i-2); elseif nx>8, tk=ttt(i-1); end
    end
    for n=1:length(xt)
       dt=mjdtocal(prm.td,xt(n)*3600);
       if mod(xt(n),tk)~=0, xl{n}='';
       elseif tk<1/60,xl{n}=sprintf('%2d:%02d:%02.0f',dt(4:6));
       elseif tk<1,   xl{n}=sprintf('%2d:%02d',       dt(4:5));
       elseif tk<24,  xl{n}=sprintf('%d/%d %2d:%02d', dt(2:5));
       elseif tk<240, xl{n}=sprintf('%d/%d',          dt(2:3));
       else xl{n}=sprintf('%02d/%d/%d',mod(dt(1),100),dt(2:3));
       end
    end
end
set(ax,'xtick',xt,'xticklabel',xl)

% draw sky map -----------------------------------------------------------------
function skymap(varargin)
back=1; mask=0; dirs=0; rot=0;
hold on, axis equal, axis([-90,90,-90,90]), axis off
n=1; va={};
while n<=nargin
    switch varargin{n}
    case 'noback', back=0; n=n+1;
    case 'mask',   mask=1; n=n+1;
    case 'dirs',   dirs=1; n=n+1;
    case 'rotate', rot=varargin{n+1}; n=n+2;
    otherwise va={va{:},varargin{n}}; n=n+1; end
end
fname=get(gca,'fontname');
fsize=get(gca,'fontsize');

if back
    patch(90*sin(0:pi/36:2*pi),90*cos(0:pi/36:2*pi),'w','linestyle','none')
end
if mask
    c=get(gca,'color');
    if strcmp(c,'none'), c=get(gcf,'color'); end
    if strcmp(c,'none'), c=get(0,'defaultuicontrolbackgroundcolor'); end
    patch(90*[sin(0:pi/36:2*pi),-1.01,-1.01,1.01,1.01],...
          90*[cos(0:pi/36:2*pi),1.01,-1.01,-1.01,1.01],c,'linestyle','none');
end
for el=0:15:90
    if el==0, s='k'; else s='k:'; end
    plot((90-el)*sin(0:pi/36:2*pi),(90-el)*cos(0:pi/36:2*pi),s)
    if el>0
        text((90-el)*sin(rot*pi/180),(90-el)*cos(rot*pi/180),[' ',num2str(el)],...
             'fontname',fname,'fontsize',fsize,va{:})
    end
end
label='NESW';
for az=0:30:330
    if mod(az,90)==0&dirs, s='k'; else s='k:'; end
    azr=az+rot;
    plot([0,90*sin(azr*pi/180)],[0,90*cos(azr*pi/180)],s)
    if dirs&mod(az,90)==0, str=label(az/90+1); else str=num2str(az); end
    text(95*sin(azr*pi/180),95*cos(azr*pi/180),str,...
         'horizontal','center','fontname',fname,'fontsize',fsize,va{:})
end

% plot in skymap ---------------------------------------------------------------
function h=skyplot(azel,varargin)
if isempty(azel), h=[]; return, end
h=plot((90-azel(:,2)*180/pi).*sin(azel(:,1)),...
       (90-azel(:,2)*180/pi).*cos(azel(:,1)),varargin{:});

% patch in skymap --------------------------------------------------------------
function h=skypatch(azel,c,varargin)
if isempty(azel), h=[]; return, end
h=patch((90-azel(:,2)*180/pi).*sin(azel(:,1)),...
        (90-azel(:,2)*180/pi).*cos(azel(:,1)),c,varargin{:});

% text in skymap ---------------------------------------------------------------
function h=skytext(azel,str,varargin)
fname=get(gca,'fontname');
fsize=get(gca,'fontsize');
h=text((90-azel(:,2)*180/pi).*sin(azel(:,1)),...
       (90-azel(:,2)*180/pi).*cos(azel(:,1)),str,'fontname',fname,...
       'fontsize',fsize,varargin{:});

% corn patch in skymap ---------------------------------------------------------
function h=skycpatch(azel,aw,ew,c,varargin)
az=(azel(1)-aw/2:pi/72:azel(1)+aw/2)';
el=azel(2)+[-ew/2;ew/2];
az0=repmat(azel(1)-aw/2,length(el),1);
az1=repmat(azel(1)+aw/2,length(el),1);
el0=repmat(azel(2)-ew/2,length(az),1);
el1=repmat(azel(2)+ew/2,length(az),1);
azel=[az,el0;az1,el;flipud(az),el1;az0,flipud(el)];
h=skypatch(azel,c,'linestyle','none',varargin{:});

% colorbar vertical ------------------------------------------------------------
function h=colorbarv(pos,range,label,varargin)
ax=gca;
h1=axes('position',pos,'ylim',range,'xtick',[],'ytick',[]);
cmap=get(gcf,'colormap');
s=(range(2)-range(1))/size(cmap,1);
for n=1:size(cmap,1)
    patch([0,1,1,0],range(1)+[n-1,n-1,n,n]*s,n,'edgecolor','none');
end
h2=axes('position',pos,'ylim',range,'box','on','xtick',[],'yaxislocation','right',...
        'color','none',varargin{:});
if nargin>=3, ylabel(label); end
axes(ax);
h=[h1;h2];

% colorbar horizontal ----------------------------------------------------------
function h=colorbarh(pos,range,label,varargin)
ax=gca;
h1=axes('position',pos,'xlim',range,'xtick',[],'ytick',[]);
cmap=get(gcf,'colormap');
s=(range(2)-range(1))/size(cmap,1);
for n=1:size(cmap,1)
    patch(range(1)+[n-1,n-1,n,n]*s,[0,1,1,0],n,'edgecolor','none');
end
h2=axes('position',pos,'xlim',range,'box','on','ytick',[],'color','none',varargin{:});
if nargin>=3, xlabel(label); end
axes(ax);
h=[h1;h2];

% plot bargraph ----------------------------------------------------------------
function h=barplot(x,labels,varargin)
[fname,fsize]=gut('getfont','g');
if size(x,1)>1, xx=x; else xx=[x;repmat(nan,1,size(x,2))]; end
if size(xx,2)>1, w=1; else w=0.5; end
h=bar(xx,w);
xl=[0.5,size(x,1)+0.5];
set(gca,'position',[0.1,0.16,0.84,0.74],'box','on','xgrid','on','ygrid','on',...
    'xlim',xl,'xtick',1.5:size(x,1)-0.5,'xticklabel',{},...
    'fontname',fname,'fontsize',fsize,varargin{:});
yl=get(gca,'ylim'); ys=yl(2)-yl(1);
for n=1:length(labels)
    text(n,yl(1)-ys*0.01,[labels{n},' '],'rotation',90,'horizontalalignment','right',...
         'fontname',fname,'fontsize',fsize);
end
set(get(gca,'xlabel'),'verticalalignment','top','position',[mean(xl),yl(1)-ys*0.13]);

% colormap dialog --------------------------------------------------------------
function [ok,cmap]=cmapdlg(tag,title,cmap)
h=gut('newdlg',tag,title,[240,90]); axis off;
if ~isempty(cmap), set(h,'colormap',cmap); end
list={'','JET','HSV','GRAY','HOT','FLAG','SPRING','SUMMER','AUTUMN','WINTER',...
      'COPPER','COOL','PINK','BONE'};
gut('newtext','',[15,30,100,20],'Color Map');
gut('newpopm','cmap',[142,34,92,20],list,lower(list),[mfilename,' cbCmap']);
colorbarh([0.04,0.68,0.92,0.2],[0,1],'','xtick',[]);
gut('newokcancelbtn','',[53,4,180,22]);
ok=gut('waitok'); if ~ok, return, end
cmap=get(h,'colormap');
close;

function cbCmap
cmap=gut('getsel','cmap');
if ~isempty(cmap), set(gcf,'colormap',feval(cmap)); end
