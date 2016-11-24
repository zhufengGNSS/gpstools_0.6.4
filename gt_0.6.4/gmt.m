function varargout=gmt(varargin)
%-------------------------------------------------------------------------------
% [system] : GpsTools
% [module] : map common libraries for GpsTools
% [func]   : map common libraries for GpsTools
% [argin]  : function,argins
%              'mmap'(,opts)          : set map properties
%                  (call before calling mcoast,mgrid,mplot,mpatch,mtext,tetlim,
%                   xytoll,lltoxy)
%                  opts : options
%                      'proj',proj    : map projection
%                          'eq'       : equidistant cylindrical
%                          'mercat'   : mercator cylindrical
%                          'miller'   : miller cylindrical
%                          'lambert'  : lambert conformal conic(phi1/2=30N/60N)
%                          'ortho'    : orthographic
%                          'polarn'   : polar stereographic north
%                          'polars'   : polar stereographic south
%                          'mollw'    : mollweide
%                      'cent',cent    : map center longitude/latitude (deg)
%                      'base',base    : map base longitude/latitude (deg)
%                      'scale',scale  : map scale
%                      'pos',pos      : map axis position
%                      'zfact',zfact  : z factor
%              'mcoast'(,opts)        : draw coast lines
%                  opts : options
%                      'lcolor',lcolor: land color
%                      'scolor',scolor: sea/lake color
%                      'ccolor',ccolor: coastline color
%                      'zpos',zpos    : z position
%              'seamask'(,opts)       : draw sea mask
%                  opts : options
%                      'scolor',mcolor: sea/lake color
%                      'zpos',zpos    : z position
%              'mgrid'(,opts)         : draw mesh grid and labels
%                  opts : options
%                      'gint',gint    : mesh grid interval (deg) (0:no grid)
%                      'lint',lint    : label interval (deg) (0:no label)
%                      'color',color  : line color
%                      'ltype',ltype  : line style
%                      'zpos',zpos    : z position
%              'mplot',lon,lat,lc,(,opts) : draw plot on map
%                  lon,lat : longitude/latitude (deg)
%                  lc      : line color
%                  opts    : options for plot
%              'mpatch',lon,lat,fc,lc,(,opts): draw patch on map
%                  lon,lat : longitude/latitude (deg)
%                  fc      : face color
%                  lc      : line color
%                  opts    : options for patch
%              'mtext',lon,lat,str,(,opts) : draw text on map
%                  lon,lat : longitude/latitude (deg)
%                  str     : text string
%                  opts    : options for text
%              'getlim'    : get x/y limit of map
%              'normlon',lon : normalize longitude in [-179.999,180]
%                  lon     : longitude (deg)
%              'lltoxy',lon,lat(,h,prm) : transform lat/lon to map cooridinate
%                  lon,lat : longitude/latitude (deg)
%                  h       : height
%                  prm     : map parameter
%                      prm.proj : map projection
%                      prm.cent : map center longitude/latitude (deg)
%                      prm.base : map base longitude/latitude (deg)
%                      prm.scale: map scale
%              'xytoll',x,y(,prm)     : transform map cooridinate to lat/lon
%                  x,y     : postion on map
%                  prm     : map parameter
%              'gridtoll',x,y,gprm    : transform grid coorinate to lat/lon
%                  x,y     : postion on grid
%                  gprm    : grid parameters struct (see readgpv.m)
%              'lltogrid',lat,lon,gprm : transform lat/lon grid coorinate
%                  lon,lat : longitude/latitude (deg)
%                  gprm    : grid parameters struct (see readgpv.m)
%              'gridtouv',x,y,ux,vx,gprm : transform grid uv to east/north uv
%                  x,y     : postion on grid
%                  u,v     : u/v-component on grid
%                  gprm    : grid parameters struct (see readgpv.m)
%              'marrow',lon,lat,vx,vy,(opts) : draw arrow (vx,vy:size(deg))
%              'setz',h,z             : set z position
%                  h       : object handles
%                  z       : z position
% [argout] :   'mcoast','mgrid','mplot','mpatch','mtext'
%                  h       : object handles
%              'getlim'
%                  xl,yl   : x/y limits of map
%              'normlon'
%                  lon     : normalized longitude (deg)
%              'lltoxy'
%                  x,y     : x/y position on map
%                  z       : z position for ortho, valid flag(>0) for lambert
%              'xytoll'
%                  lon,lat : longitude/latitude (deg)
%              'gridtoll'
%                  lon,lat : longitude/latitude (deg)
%              'lltogrid'
%                  x,y     : x/y position on grid
%              'gridtouv'
%                  u,v     : east/north compornent on map
% [note]   :
% [version]: $Revision: 12 $ $Date: 2008-11-25 10:02:15 +0900 (火, 25 11 2008) $
%            Copyright(c) 2004-2006 by T.Takasu, all rights reserved
% [history]: 04/10/29  0.1  new
%            05/03/14  0.2  add mmap
%            05/05/03  0.3  add lambert projection
%            05/05/08  0.4  add function gridtoll,lltogrid,gridtouv,setz
%                           add clipping of coastlines
%            05/09/29  0.5  add mollweide projection
%            05/11/23  0.6  add function seamask
%-------------------------------------------------------------------------------
if nargout==0, feval(varargin{:});
else [varargout{1:nargout}]=feval(varargin{:}); end

% set map properties -----------------------------------------------------------
function mmap(varargin)
global g_gmt
g_gmt.proj='eq'; g_gmt.cent=[0,0]; g_gmt.base=[0,0]; g_gmt.scale=1;
pos=[0.06,0.09,0.88,0.82]; zfact=0.1; n=1; va={};
while n<=nargin
    switch varargin{n}
    case 'proj', g_gmt.proj =varargin{n+1}; n=n+2; % map projection
    case 'cent', g_gmt.cent =varargin{n+1}; n=n+2; % center lon/lat (deg)
    case 'base', g_gmt.base =varargin{n+1}; n=n+2; % base lon/lat (deg)
    case 'scale',g_gmt.scale=varargin{n+1}; n=n+2; % scale
    case 'pos',  pos=varargin{n+1}; n=n+2;         % axes position
    case 'zfact',zfact=varargin{n+1}; n=n+2;       % z factor
    otherwise, va={va{:},varargin{n:n+1}}; n=n+2; end
end
hold on
set(gca,'units','normalized');
set(gca,'position',pos);
set(gca,'units','pixel');
xl=[-1,1]; yl=[-1,1]; p=get(gca,'position'); r=(p(4)-p(2))/(p(3)-p(1));
switch g_gmt.proj
case 'ortho', xl=xl/r; zl=1;
case {'lambert','polarn','polars'}, xl=xl/r; zl=zfact;
otherwise, yl=yl*r; zl=zfact;
end
set(gca,'units','normalized','xlim',xl,'ylim',yl,'xtick',[],'ytick',[],...
    'plotboxaspectratio',[xl(2),yl(2),zl],va{:});

% get coast lines --------------------------------------------------------------
function [lon,lat,type]=getcoast % (type=1:lake,0:land)
global g_gmt, if isempty(g_gmt), error('set map properties first'); end
persistent pdb p info
switch g_gmt.proj
case 'ortho'
    if     g_gmt.scale>8,  db='gshhs_i.mat';
    elseif g_gmt.scale>2,  db='gshhs_l.mat';
    else                   db='gshhs_w.mat'; end
case {'lambert','polarn','polars'}
    if     g_gmt.scale>6,  db='gshhs_i.mat';
    elseif g_gmt.scale>1.5,db='gshhs_l.mat';
    else                   db='gshhs_w.mat'; end
otherwise
    if     g_gmt.scale>16, db='gshhs_i.mat';
    elseif g_gmt.scale>4,  db='gshhs_l.mat';
    else                   db='gshhs_w.mat'; end
end
if isempty(pdb)|~strcmp(db,pdb)
    [dirs,f]=fileparts(which(mfilename));
    file=fullfile(dirs,'data',db);
    load(file); info=[p.inf]'; pdb=db;
    info(info(:,5)==360,5)=359.999;
end
lon={}; lat={}; type=[]; n=0; 
for i=clipmap(info(:,4:5),info(:,6:7))';
    n=n+1; lon{n}=p(i).lon; lat{n}=p(i).lat;
    if mod(info(i,3),2)==0, type(n)=0;
    else type(n)=1; lon{n}=flipud(lon{n}); lat{n}=flipud(lat{n}); end
end

% draw coast lines -------------------------------------------------------------
function h=mcoast(varargin)
global g_gmt, if isempty(g_gmt), error('set map properties first'); end
h=[]; ccolor='k'; lcolor='none'; scolor='none'; zpos=[]; n=1; va={};
while n<=nargin
    switch varargin{n}
    case 'lcolor', lcolor=varargin{n+1}; n=n+2; % land color
    case 'scolor', scolor=varargin{n+1}; n=n+2; % sea color
    case 'ccolor', ccolor=varargin{n+1}; n=n+2; % coastline color
    case 'zpos',   zpos  =varargin{n+1}; n=n+2; % z position
    otherwise, va={va{:},varargin{n:n+1}}; n=n+2; end
end
if strcmp(lcolor,'none'), scolor='none';
elseif strcmp(scolor,'none')
    if strcmp(get(gca,'visible'),'on'), scolor=get(gca,'color');
    else scolor=get(gcf,'color'); end
end
switch g_gmt.proj % map edges
case 'ortho'
    t=(0:5:360)*pi/180; x=cos(t)*g_gmt.scale; y=sin(t)*g_gmt.scale;
case 'mollw'
    t=(0:5:360)*pi/180; x=cos(t)*g_gmt.scale; y=sin(t)*g_gmt.scale/2;
otherwise
    [xl,yl]=getlim; x=xl([1,1,2,2,1]); y=yl([1,2,2,1,1]);
end
if ~strcmp(scolor,'none')
    h=[h;patch(x,y,scolor,'edgecolor','none');plot(x,y,'k')];
end
[lon,lat,type]=getcoast;
for n=1:length(lon)
    if type(n), fc=lcolor; else fc=scolor; end
    h=[h;mpatch(lon{n},lat{n},fc,ccolor,va{:})];
end
if ~strcmp(ccolor,'none'), h=[h;plot(x,y,'color',ccolor)]; end
if ~isempty(zpos), setz(h,zpos); end

% clip map extent --------------------------------------------------------------
function i=clipmap(lone,late)
global g_gmt, if isempty(g_gmt), error('set map properties first'); end
[xl,yl]=getlim; [lon,lat]=xytoll(xl([1,1,2,2]),yl([1,2,1,2]));
switch g_gmt.proj
case {'eq','mercat','miller'}
    lat=[minn(lat),maxn(lat)];
    lon=[g_gmt.cent(1),180/g_gmt.scale];
case 'ortho'
    if any(isnan(lon)|isnan(lat))
        lat=[-90,90]; lon=[0,180];
    else
        p=asin(1/g_gmt.scale)*180/pi;
        lat=[minn([lat;g_gmt.cent(2)-p]),maxn([lat;g_gmt.cent(2)+p])];
        lon=[g_gmt.cent(1),maxn(abs(normlon(lon-g_gmt.cent(1))))];
    end
case 'lambert'
    [x,y]=lltoxy(g_gmt.base(1),90);
    [lonn,latn]=xytoll(x,yl(2)); if isnan(latn), latn=90; end
    lats=minn(lat); if isempty(lats), lats=-90; end
    lat=[lats,latn]; lon=[g_gmt.cent(1),maxn(abs(normlon(lon-g_gmt.cent(1))))];
case 'polarn'
    lat=[minn(lat),90]; lon=[-180,180];
case 'polars'
    lat=[-90,maxn(lat)]; lon=[-180,180];
otherwise
    lat=[-90,90]; lon=[-180,180];
end
if length(lon)<2|length(lat)<2, i=[]; return, end
d=normlon(lone-lon(1));
i=find((any(abs(d(:,1:2))<=lon(2),2)|d(:,2)<d(:,1)|(d(:,1)<0&0<d(:,2)))&...
       lat(1)<=late(:,2)&late(:,1)<=lat(2));

% max/min without nan ----------------------------------------------------------
function x=maxn(x), x=max(x(~isnan(x)));
function x=minn(x), x=min(x(~isnan(x)));

% draw seamask -----------------------------------------------------------------
function h=seamask(varargin)
global g_gmt, if isempty(g_gmt), error('set map properties first'); end
h=[]; scolor='w'; zpos=[]; n=1; va={};
while n<=nargin
    switch varargin{n}
    case 'scolor', scolor=varargin{n+1}; n=n+2; % sea/lake color
    case 'zpos',   zpos  =varargin{n+1}; n=n+2; % z position
    otherwise, va={va{:},varargin{n:n+1}}; n=n+2; end
end
[lon,lat,type]=getcoast;
[xl,yl]=getlim; xs={}; ys={}; xa=xl; ya=yl;
for n=1:length(lon)
    if ~type(n)
        h=[h;mpatch(lon{n},lat{n},scolor,'none',va{:})]; % lake
    else
        [x,y,z]=lltoxy(lon{n},lat{n});
        [x,y]=sepseg(x,y,z,2);
        for m=1:length(x)
            a=[min(x{m}),max(x{m}),min(y{m}),max(y{m})];
            if ~isempty(a)&xl(1)<=a(2)&a(1)<=xl(2)&yl(1)<=a(4)&a(3)<=yl(2)
                xs={xs{:},x{m}}; ys={ys{:},y{m}};
                xa=[min([xa(1),a(1)]),max([xa(2),a(2)])];
                ya=[min([ya(1),a(3)]),max([ya(2),a(4)])];
            end
        end
    end
end
xs={xs{:},[xa(1);xa(2)+0.001]}; ys={ys{:},[ya(1);ya(1)]}; % add bottom line
for n=1:length(ys)-1
    [yy,i]=min(ys{n}); % bottom point of map
    mm=[]; jj=[]; pp=[];
    for m=n+1:length(ys)
        [j,p]=crossp(xs{m},ys{m},xs{n}(i),ys{n}(i));
        if ~isempty(j), mm=[mm;m]; jj=[jj;j]; pp=[pp;p]; end
    end
    if ~isempty(mm)
        [yy,k]=max(pp(:,2)); m=mm(k);
        [xs{m},ys{m}]=catmap(xs{n},ys{n},xs{m},ys{m},i,jj(k),pp(k,:));
    end
end
xs=[xs{end};xa([2,1,1])'];
ys=[ys{end};ya([2,2,1])'];
h=[h;patch(xs,ys,scolor,'edgecolor','none',va{:})];
if ~isempty(zpos), setz(h,zpos); end

% crossing point with lower vertical half-line ---------------------------------
function [i,p]=crossp(x,y,xx,yy)
i=[]; p=[]; if xx<min(x)|max(x)<xx, return, end
for j=find(x(1:end-1)<=xx&xx<x(2:end))'
    i=[i;j]; a=(xx-x(j))/(x(j+1)-x(j));
    p=[p;x(j)*(1-a)+x(j+1)*a,y(j)*(1-a)+y(j+1)*a];
end
if isempty(i), return, end
j=find(p(:,2)<=yy); [ymax,k]=max(p(j,2));
i=i(j(k)); p=p(j(k),:);

% connect maps -----------------------------------------------------------------
function [x,y]=catmap(x1,y1,x2,y2,i,j,p)
x=[x2(1:j);p(1);x1(i:end);x1(1:i);p(1);x2(j+1:end)];
y=[y2(1:j);p(2);y1(i:end);y1(1:i);p(2);y2(j+1:end)];

% draw mesh grid and labels ----------------------------------------------------
function h=mgrid(varargin)
global g_gmt, if isempty(g_gmt), error('set map properties first'); end
h=[]; gint=15; lint=30; color='k'; ltype=':'; zpos=[]; n=1; va={};
while n<=nargin
    switch varargin{n}
    case 'gint',  gint =varargin{n+1}; n=n+2; % mesh grid interval (deg)
    case 'lint',  lint =varargin{n+1}; n=n+2; % lat/lon label interval (deg)
    case 'color', color=varargin{n+1}; n=n+2; % mesh grid color
    case 'ltype', ltype=varargin{n+1}; n=n+2; % line type
    case 'zpos',  zpos =varargin{n+1}; n=n+2; % z position
    otherwise, va={va{:},varargin{n:n+1}}; n=n+2; end
end
[xl,yl]=getlim;
if g_gmt.scale<3&~strcmp(g_gmt.proj,'mollw'), eint=5; else eint=2; end
if gint>0&~strcmp(color,'none') % mesh grid
    for lon=-180:gint:180-gint
        switch g_gmt.proj
        case {'eq','miller'}
            lats=[-89.999,89.999];
        case 'mercat'
            lats=[-85,85];
        case {'ortho','lambert','polarn','polars'}
            if mod(lon,90)==0, lats=-90:eint:90; else lats=gint-90:eint:90-gint; end
        otherwise
            lats=-90:eint:90;
        end
        h=[h;mplot(repmat(lon,size(lats)),lats,color,'linestyle',ltype,va{:})];
    end
    for lat=-90:gint:90
        if lat==0, s='-'; else s=ltype; end
        switch g_gmt.proj
        case {'eq','mercat','miller'}
            [x,y]=lltoxy(0,lat); h=[h;plot(xl,[y,y],s,'color',color,va{:})];
        otherwise
            lons=-180:eint:180;
            h=[h;mplot(lons,repmat(lat,size(lons)),color,'linestyle',s,va{:})];
        end
    end
end
if lint>0 % lat/lon labels
    xs={'W','','E'}; ys={'S','','N'};
    for lon=-180+lint:lint:180
        s=[num2str(abs(lon)),'\circ',xs{sign(lon)+2}];
        switch g_gmt.proj
        case {'eq','mercat','miller'},
            [x,y]=lltoxy(lon,0);
            if xl(1)<=x&x<=xl(2), h=[h;PlotText(x,yl(1),s,2,3)]; end
        case 'lambert'
            lat=-60:g_gmt.cent(2);
            [x,y,z]=lltoxy(repmat(lon,size(lat)),lat);
            i=find(xl(1)<x&x<xl(2)&yl(1)<y&y<yl(2)&z>0);
            [yy,j]=min(y(i)); i=i(j);
            if ~isempty(i)&yy<yl(1)*0.8, h=[h;PlotText(x(i),yl(1),s,2,3)]; end
        otherwise
            [x,y,z]=lltoxy(lon,0.001);
            if xl(1)<=x&x<=xl(2)&z>=0, h=[h;PlotText(x,y,s,2,2)]; end
        end
    end
    for lat=-90+lint:lint:90-lint
        s=[num2str(abs(lat)),'\circ',ys{sign(lat)+2}];
        switch g_gmt.proj
        case {'eq','mercat','miller'}
            [x,y,z]=lltoxy(0,lat);
            if yl(1)<=y&y<=yl(2), h=[h;PlotText(xl(1),y,[' ',s],1,2)]; end
        otherwise
            [x,y,z]=lltoxy(g_gmt.cent(1),lat);
            if yl(1)<=y&y<=yl(2)&z>=0, h=[h;PlotText(x,y,s,2,2)]; end
        end
    end
end
if ~isempty(zpos), setz(h,zpos); end

% plot text --------------------------------------------------------------------
function h=PlotText(x,y,s,hp,vp)
fn=get(gca,'fontname'); fs=get(gca,'fontsize'); 
hpos={'left','center','right'}; vpos={'top','middle','bottom'};
h=text(x,y,s,'horizontalalignment',hpos{hp},'verticalalignment',vpos{vp},...
       'fontname',fn,'fontsize',fs);

% plot in map ------------------------------------------------------------------
function h=mplot(lon,lat,lc,varargin)
h=[]; if strcmp(lc,'none'), return, end
[x,y,z]=lltoxy(lon,lat);
[xs,ys]=sepseg(x,y,z,1);
for n=1:length(xs), h=[h;plot(xs{n},ys{n},'color',lc,varargin{:})]; end

% patch in map -----------------------------------------------------------------
function h=mpatch(lon,lat,fc,lc,varargin)
h=[];
[x,y,z]=lltoxy(lon,lat);
[xs,ys]=sepseg(x,y,z,2);
if strcmp(fc,'none'), fc=nan; end
for n=1:length(xs), h=[h;patch(xs{n},ys{n},fc,'edgecolor',lc,varargin{:})]; end

% text in map ------------------------------------------------------------------
function h=mtext(lon,lat,str,varargin)
h=[];
[x,y,z]=lltoxy(lon,lat);
xl=xlim; yl=ylim;
if z>=0&xl(1)<=x&x<=xl(2)&yl(1)<=y&y<=yl(2), h=text(x,y,str,varargin{:}); end

% set z position ---------------------------------------------------------------
function setz(h,z)
if isempty(z), return, end
for n=1:length(h)
    switch get(h(n),'type')
    case {'line','patch'}, set(h(n),'zdata',repmat(z,size(get(h(n),'xdata'))));
    case 'text', p=get(h(n),'position'); set(h(n),'position',[p(1:2),z]);
    end
end

% get x/y limit ----------------------------------------------------------------
function [xl,yl]=getlim
global g_gmt, if isempty(g_gmt), error('set map properties first'); end
xl=get(gca,'xlim'); yl=get(gca,'ylim');
switch g_gmt.proj
case {'eq','mercat','miller'}
    if strcmp(g_gmt.proj,'mercat'), lat=85; else lat=90; end
    [x,y]=gmt('lltoxy',[0,0],[-lat,lat]);
    xl=[max(-g_gmt.scale,xl(1)),min(g_gmt.scale,xl(2))];
    yl=[max(y(1),yl(1)),min(y(2),yl(2))];
case 'lambert'
    [x,y]=lltoxy(0,90); yl(2)=min(yl(2),y);
end

% transform lat/lon to map coordinate ------------------------------------------
function [x,y,z]=lltoxy(lon,lat,h,prm)
global g_gmt
if nargin<3, h=zeros(size(lon)); end
if nargin<4
    if isempty(g_gmt), error('set map properties first'); end
    prm=g_gmt;
end
lat=lat*pi/180;
switch prm.proj
case {'eq','mercat','miller'}
    lat0=prm.cent(2)*pi/180;
    x=normlon(lon-prm.cent(1))/180; 
    switch prm.proj
    case 'eq',     y=(lat-lat0)/pi;
    case 'mercat', y=(atanh(sin(lat))-atanh(sin(lat0)))/pi;
    case 'miller', y=(atanh(sin(lat*0.8))-atanh(sin(lat0*0.8)))/0.8/pi;
    end
    z=h;
case 'ortho'
    lat0=prm.cent(2)*pi/180;
    lon=(lon-prm.cent(1))*pi/180;
    sinl=sin(lon); cosl=cos(lon); sinp=sin(lat); cosp=cos(lat); 
    x=(1+h).*cosp.*sinl;
    y=(1+h).*(cos(lat0)*sinp-sin(lat0)*cosp.*cosl);
    z=(1+h).*(sin(lat0)*sinp+cos(lat0)*cosp.*cosl);
case 'lambert'
    lat0=min(max(prm.cent(2),-60),89.999)*pi/180;
    lon0=normlon(prm.cent(1)-prm.base(1))*pi/180;
    lon=normlon(lon-prm.base(1))*pi/180;
    x=repmat(nan,size(lon)); y=x; z=x;
    i=find(-60*pi/180<lat&lat<=pi/2);
    n=0.7155668; f=1.7930256; % phi1=30deg,phi2=60deg
    p0=f./tan(pi/4+lat0/2).^n;
    p=f./tan(pi/4+lat(i)/2).^n;
    z(i)=p.*cos(n*lon(i));
    x(i)=p.*sin(n*lon(i))-p0*sin(n*lon0);
    y(i)=p0*cos(n*lon0)-z(i);
case 'polarn'
    lon0=(prm.cent(1)-prm.base(1))*pi/180;
    scale=(1+sin(60*pi/180))*prm.scale;
    r0=scale*tan(pi/4-prm.cent(2)*pi/360);
    x=repmat(nan,size(lon)); y=x; z=x;
    i=find(lat>-70*pi/180);
    r=scale*tan(pi/4-lat(i)/2);
    lon=normlon(lon(i)-prm.base(1))*pi/180;
    x(i)=r.*sin(lon)-r0*sin(lon0);
    y(i)=r0*cos(lon0)-r.*cos(lon);
    z(i)=h(i);
case 'polars'
    lon0=(prm.cent(1)-prm.base(1))*pi/180;
    scale=(1+sin(60*pi/180))*prm.scale;
    r0=scale*tan(pi/4+prm.cent(2)*pi/360);
    x=repmat(nan,size(lon)); y=x; z=x;
    i=find(lat<80*pi/180);
    r=scale*tan(pi/4+lat(i)/2);
    lon=normlon(lon(i)-prm.base(1))*pi/180;
    x(i)=r.*sin(lon)-r0*sin(lon0);
    y(i)=r.*cos(lon)-r0*cos(lon0);
    z(i)=h(i);
case 'mollw'
    x=normlon(lon-prm.cent(1))/180.*sqrt(1-(2*lat/pi).^2);
    y=lat/pi;
    z=h;
end
x=x*prm.scale; y=y*prm.scale; z=z*prm.scale;

% transform map coordinate to lat/lon ------------------------------------------
function [lon,lat]=xytoll(x,y,prm)
global g_gmt
if nargin<3
    if isempty(g_gmt), error('set map properties first'); end
    prm=g_gmt;
end
x=x(:)/prm.scale; y=y(:)/prm.scale;
lon=repmat(nan,length(x),1); lat=lon;
switch prm.proj
case {'eq','mercat','miller'}
    y=y*pi; y0=prm.cent(2)*pi/180;
    switch prm.proj
    case 'mercat', y=2*(atan(exp((y+y0)))-atan(exp(y0)))+y0;
    case 'miller', y=2*(atan(exp((y+y0)*0.8))-atan(exp(y0*0.8)))/0.8+y0;
    otherwise, y=y+y0; end
    lon=x*180+prm.cent(1);
    lat=y*180/pi;
case 'ortho'
    cosp=cos(prm.cent(2)*pi/180); sinp=sin(prm.cent(2)*pi/180);
    r=x.^2+y.^2;
    i=find(r<=1);
    z=sqrt(1-r(i));
    lon(i)=atan2(x(i),cosp*z-sinp*y(i))*180/pi+prm.cent(1);
    lat(i)=asin(cosp*y(i)+sinp*z)*180/pi;
case 'lambert'
    lon0=normlon(prm.cent(1)-prm.base(1))*pi/180; 
    lat0=min(max(prm.cent(2),-60),89.999)*pi/180;
    n=0.7155668; f=1.7930256; % phi1=30deg,phi2=60deg
    p0=f./tan(pi/4+lat0/2).^n;
    x=x+p0*sin(n*lon0);
    y=p0*cos(n*lon0)-y;
    i=find(y>=0);
    lat(i)=atan((f./sqrt(x(i).^2+y(i).^2)).^(1/n))*360/pi-90;
    lon(i)=prm.base(1)+atan(x(i)./y(i))/n*180/pi;
case 'polarn'
    lon0=(prm.cent(1)-prm.base(1))*pi/180;
    scale=(1+sin(60*pi/180))*prm.scale;
    r0=scale*tan(pi/4-prm.cent(2)*pi/360);
    x=x+r0*sin(lon0);
    y=r0*cos(lon0)-y;
    lon=normlon(atan2(x,y)*180/pi+prm.base(1));
    lat=90-atan(sqrt(x.^2+y.^2)/scale)*360/pi;
case 'polars'
    lon0=(prm.cent(1)-prm.base(1))*pi/180;
    scale=(1+sin(60*pi/180))*prm.scale;
    r0=scale*tan(pi/4+prm.cent(2)*pi/360);
    x=x+r0*sin(lon0);
    y=y+r0*cos(lon0);
    lon=normlon(atan2(x,y)*180/pi+prm.base(1));
    lat=atan(sqrt(x.^2+y.^2)/scale)*360/pi-90;
case 'mollw'
    i=find(x.^2+(2*y).^2<=1);
    lat(i)=180*y(i);
    lon(i)=prm.cent(1)+180*x(i)./sqrt(1-(lat(i)/90).^2);
end
lon=normlon(lon);

% separate segments (opt=1:for plot,opt=2:for patch) --------------------------
function [xs,ys]=sepseg(x,y,z,opt)
global g_gmt, if isempty(g_gmt), error('set map properties first'); end
x=x(:); y=y(:); z=z(:); xs{1}=[]; ys{1}=[]; 
switch g_gmt.proj
case {'eq','mercat','miller'}
    v=abs(diff(x))>g_gmt.scale; if all(~v), xs{1}=x; ys{1}=y; return, end
    i=1; n=1; xn=[]; yn=[];
    for j=find(v)'
        [xa,ya,xb,yb]=AddSideP(x(j:j+1),y(j:j+1),g_gmt.scale);
        xs{n}=[xn;x(i:j);xa]; ys{n}=[yn;y(i:j);ya];
        xn=xb; yn=yb; i=j+1; n=n+1;
    end
    if opt==2
        xs{1}=[xn;x(i:end);xs{1}]; ys{1}=[yn;y(i:end);ys{1}];
        [xs,ys]=AddSideE(xs,ys);
    else
        xs{n}=[xn;x(i:end)]; ys{n}=[yn;y(i:end)];
    end
case 'ortho'
    v=z>=0; if all(v), xs{1}=x; ys{1}=y; return, elseif all(~v), return, end
    i=1; n=1; xn=[]; yn=[];
    for j=find(diff(v))'
        [xa,ya]=AddArcP(x(j:j+1),y(j:j+1),z(j:j+1),g_gmt.scale);
        if z(j)>=0, xs{n}=[xn;x(i:j);xa]; ys{n}=[yn;y(i:j);ya]; n=n+1; i=[];
        else xn=xa; yn=ya; i=j+1; end
    end
    if opt==2
        if ~isempty(i), xs{1}=[xn;x(i:end);xs{1}]; ys{1}=[yn;y(i:end);ys{1}]; end
        [xs,ys]=AddArcE(xs,ys,g_gmt.scale,g_gmt.scale);
    else
        if ~isempty(i), xs{n}=[xn;x(i:end)]; ys{n}=[yn;y(i:end)]; end
    end
case 'lambert'
    v=z>=0; if all(v), xs{1}=x; ys{1}=y; return, elseif all(~v), return, end
    i=1; n=1; xn=[]; yn=[];
    for j=find(diff(v))'
        [xa,ya]=AddUpperP(x(j:j+1),y(j:j+1),z(j:j+1));
        if z(j)>=0, xs{n}=[xn;x(i:j);xa]; ys{n}=[yn;y(i:j);ya]; n=n+1; i=[];
        else xn=xa; yn=ya; i=j+1; end
    end
    if opt==2
        if ~isempty(i), xs{1}=[xn;x(i:end);xs{1}]; ys{1}=[yn;y(i:end);ys{1}]; end
        [xs,ys]=AddUpperE(xs,ys);
    else
        if ~isempty(i), xs{n}=[xn;x(i:end)]; ys{n}=[yn;y(i:end)]; end
    end
case 'mollw'
    xx=x./sqrt(1.0001-(2*y/g_gmt.scale).^2);
    v=abs(diff(xx))>g_gmt.scale; if all(~v), xs{1}=x; ys{1}=y; return, end
    i=1; n=1; xn=[]; yn=[];
    for j=find(v)'
        [xa,ya,xb,yb]=AddSideP(xx(j:j+1),y(j:j+1),g_gmt.scale);
        xa=xa*sqrt(1-(2*ya/g_gmt.scale)^2);
        xb=xb*sqrt(1-(2*yb/g_gmt.scale)^2);
        xs{n}=[xn;x(i:j);xa]; ys{n}=[yn;y(i:j);ya];
        xn=xb; yn=yb; i=j+1; n=n+1;
    end
    if opt==2
        xs{1}=[xn;x(i:end);xs{1}]; ys{1}=[yn;y(i:end);ys{1}];
        [xs,ys]=AddArcE(xs,ys,g_gmt.scale,g_gmt.scale/2);
    else
        xs{n}=[xn;x(i:end)]; ys{n}=[yn;y(i:end)];
    end
otherwise
    xs{1}=x; ys{1}=y;
end

% add complemental side edge points --------------------------------------------
function [xa,ya,xb,yb]=AddSideP(x,y,scale)
if x(1)>=0, a=(scale-x(1))/(x(2)+2*scale-x(1)); xa=scale; xb=-scale; 
else a=(scale+x(1))/(x(1)+2*scale-x(2)); xa=-scale; xb=scale; end
ya=y(1)*(1-a)+y(2)*a; yb=ya;

% add complemental upper edge points -------------------------------------------
function [xa,ya]=AddUpperP(x,y,z)
xa=(x(1)*z(2)-x(2)*z(1))/(z(2)-z(1));
ya=(y(1)*z(2)-y(2)*z(1))/(z(2)-z(1));

% add complemental side edge ---------------------------------------------------
function [xc,yc]=AddSideE(xs,ys)
[xl,yl]=getlim; 
for n=1:min(2,length(xs))
    xc{n}=xs{n}; yc{n}=ys{n};
    for m=n+2:2:length(xs)
        [xa,ya]=cedge(xc{n}(end),yc{n}(end),xs{m}(1),ys{m}(1),yl);
        xc{n}=[xc{n};xa;xs{m}];
        yc{n}=[yc{n};ya;ys{m}];
    end
    [xa,ya]=cedge(xc{n}(end),yc{n}(end),xc{n}(1),yc{n}(1),yl);
    xc{n}=[xc{n};xa]; yc{n}=[yc{n};ya];
end

function [xa,ya]=cedge(x1,y1,x2,y2,yl)
if x1==x2, xa=[]; ya=[];
else xa=[x1;x2]; if x1>x2, ya=[yl(1);yl(1)]; else ya=[yl(2);yl(2)]; end, end

% add complemental upper edge --------------------------------------------------
function [xc,yc]=AddUpperE(xs,ys)
[x,y]=lltoxy(0,90); xc{1}=[]; yc{1}=[];
for n=1:length(xs), xc{1}=[xc{1};x;xs{n}]; yc{1}=[yc{1};y;ys{n}]; end

% add complemental arc point for circular projection ---------------------------
function [xa,ya]=AddArcP(x,y,z,scale)
a=z(2)/(z(2)-z(1));
t=atan2(y(1)*a+y(2)*(1-a),x(1)*a+x(2)*(1-a));
xa=scale*cos(t);
ya=scale*sin(t);

% add complemental arc edges for circular projection ---------------------------
function [xc,yc]=AddArcE(xs,ys,scalex,scaley)
for n=1:length(xs)
    ts(n)=atan2(ys{n}(1)/scaley,xs{n}(1)/scalex);
    te(n)=atan2(ys{n}(end)/scaley,xs{n}(end)/scalex);
end
[ts,i]=sort(ts); te=te(i); xs=xs(i); ys=ys(i);
n=1; xc={}; yc={}; d=5*pi/180;
while any(~isnan(ts))
    xc{n}=[]; yc{n}=[]; i=min(find(~isnan(ts))); j=i;
    while 1
        k=max(find(ts<=te(j))); if isempty(k), k=max(find(~isnan(ts))); end
        if ts(k)<=te(j), t=te(j):-d:ts(k); else t=te(j)+2*pi:-d:ts(k); end
        xc{n}=[xc{n};xs{j};scalex*cos(t')];
        yc{n}=[yc{n};ys{j};scaley*sin(t')];
        ts(k)=nan; if k~=i, j=k; else break, end
    end
    n=n+1;
end

% normalize longitude in [-179.999,180] ----------------------------------------
function lon=normlon(lon), lon=lon+floor((180-lon)/360)*360;

% transform grid coordinate to lon/lat -----------------------------------------
function [lon,lat]=gridtoll(x,y,gprm)
switch gprm.type
case 0 % lat/lon grid
    lon=gprm.lon1+(x-1)*gprm.dx;
    lat=gprm.lat1-(y-1)*gprm.dy;
case 3 % lambert(JMA MSM/RSM lat1=30N,lat2=60N)
    x=(x-gprm.x0)*gprm.dx;
    y=(y-gprm.y0)*gprm.dy+7.7106E6;
    lon=1.3975*atan(x./y)*180/pi+gprm.lov;
    lat=90-360*atan(1.37003E-10*(x.^2+y.^2).^0.69875)/pi;
end

% transform lat/lon to grid coordnate ------------------------------------------
function [x,y]=lltogrid(lon,lat,gprm)
switch gprm.type
case 0 % lat/lon grid
    x=(lon-gprm.lon1)/gprm.dx+1;
    y=(gprm.lat1-lat)/gprm.dy+1;
case 3 % lambert(JMA MSM/RSM lat1=30N,lat2=60N)
    x=1.14234E7*tan(pi/4-lat*pi/360).^0.71557.*sin(0.71557*(lon-gprm.lov)*pi/180);
    y=1.14234E7*tan(pi/4-lat*pi/360).^0.71557.*cos(0.71557*(lon-gprm.lov)*pi/180);
    x=x/gprm.dx+gprm.x0;
    y=(y-7.71061E6)/gprm.dy+gprm.y0;
end

% transform grid uv to lon/lat uv ----------------------------------------------
function [u,v]=gridtouv(x,y,ux,vx,gprm)
switch gprm.type
case 0 % lat/lon grid
    u=ux; v=vx;
case 3 % lambert(JMA MSM/RSM lat1=30N,lat2=60N)
    x=x-gprm.x0; y=y-gprm.y0+7.71061E6/gprm.dy;
    d=sqrt(x.^2+y.^2);
    fs=x./d; fc=y./d; u=ux.*fc+vx.*fs; v=vx.*fc-ux.*fs;
end

% draw arrow -------------------------------------------------------------------
function h=marrow(lon,lat,vx,vy,varargin)
if isempty(lon)|isempty(lat), h=[]; return, end
[x,y]=lltoxy(lon,lat);
xl=get(gca,'xlim');
f=(xl(2)-xl(1))/100; % size factor
%s=0.4*f; % head size
s=0.2*f; % head size
r=[vx,vy]*f;
rr=norm(r); if rr<s*2, return, end
h=plot([x,x+r(1)],[y,y+r(2)],varargin{:});
g=[r(1),-r(2);r(2),r(1)]/rr*[-3,0,-3;1,0,-1]*s; % arrow head
h=[h;plot(x+r(1)+g(1,:),y+r(2)+g(2,:),'-','color',get(h,'color'))];
