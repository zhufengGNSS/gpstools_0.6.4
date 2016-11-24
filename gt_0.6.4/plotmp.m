function h=plotmp(varargin)
%-------------------------------------------------------------------------------
% [system] : GpsTools
% [module] : plot multipath profile
% [func]   : plot multipath profile
% [argin]  : file = multipath profile files
% [argout] : h    = figure handle
% [note]   :
% [version]: $Revision: 16 $ $Date: 2008-12-12 15:49:30 +0900 (金, 12 12 2008) $
%            Copyright(c) 2004-2008 by T.Takasu, all rights reserved
% [history]: 04/11/19  0.1  new
%            04/12/14  0.2  support gui interface
%            08/12/02  0.3  support uigetfile bug of matlab 7.4
%-------------------------------------------------------------------------------
if nargin>0&strncmp(varargin{1},'cb',2), feval(varargin{:})
else h=PlotMain(varargin{:}); end

% plot main --------------------------------------------------------------------
function h=PlotMain(varargin)
h=gut('newfigg','','Multipath Profile',[600,400,0,-68]);
if isempty(h), return, end
gut('newmenu','data','&Data ',{'&Read File...'},{},[mfilename,' cbRead']);
gut('newmenu','plot','&Plot',...
    {'Antenna PCV/Phase Multipath (IF)','Code Multipath (L1)','Code Multipath (L2)'},...
    {'pmp','cmp1','cmp2'},[mfilename,' cbType']);
gut('newmenu','opt','&Option',...
    {'Az/El Graph / Sky Plot','-','Show Satellite Path','Show Contour Line',...
     'Value Range...'},{},[mfilename,' cbOption']);

prm.type='pmp'; prm.gtype=0; prm.nopath=1; prm.range=[]; prm.cont=0;
data.prm=prm; set(h,'userdata',data);
if nargin>0, ReadData(varargin{1}); PlotData; end

% callback on menu read --------------------------------------------------------
function cbRead
persistent file, if isempty(file), file='*.mat'; end
[f,path]=uigetfile(file,'Read Multipath Profile');
if f<=0, return; else file=[path,f]; end
ReadData(file); PlotData;

% callback on type change ------------------------------------------------------
function cbType
data=get(gcf,'userdata');
data.prm.type=get(gcbo,'userdata');
gut('unchkmenu','plot');
gut('chkmenu','plot',data.prm.type);
set(gcf,'userdata',data);
PlotData;

% callback on option change ----------------------------------------------------
function cbOption
data=get(gcf,'userdata');
switch get(gcbo,'userdata')
case 1, if gut('togglechk',gcbo), data.prm.gtype =1; else data.prm.gtype =0; end
case 2, if gut('togglechk',gcbo), data.prm.nopath=0; else data.prm.nopath=1; end
case 3, if gut('togglechk',gcbo), data.prm.cont  =1; else data.prm.cont  =0; end
case 4
    [range,ok]=gut('rangedlg','Value Range (m)',[190,60],data.prm.range);
    if ~ok, return, end
    if range(1)<range(2), data.prm.range=range; else data.prm.range=[]; end
end
set(gcf,'userdata',data);
PlotData

% plot multipath file ----------------------------------------------------------
function PlotData
data=get(gcf,'userdata');
ep1=mjdtocal(data.td+data.time(1)/86400);
ep2=mjdtocal(data.td+data.time(end)/86400);
ti=[tstr(ep1),'-',tstr(ep2)];
switch data.prm.type
case 'cmp1', ti=[data.rcv,' Code Multipath (L1) (m) : ',ti];
case 'cmp2', ti=[data.rcv,' Code Multipath (L2) (m) : ',ti];
case 'pmp',  ti=[data.rcv,' Antenna PCV/Phase Multipath (IF) (m) : ',ti];
end
if ~isempty(data.prm.range), range=data.prm.range;
elseif strcmp(data.prm.type,'pmp'), range=[-0.05,0.05];
else range=[-0.5,0.5]; end

if data.prm.nopath, azels=[]; else azels=data.azels; end
if data.prm.gtype
    switch data.prm.type
    case 'cmp1', PlotMpMapA(data.mpcc1,data.mpcs1,azels,range,data.prm.cont,ti)
    case 'cmp2', PlotMpMapA(data.mpcc2,data.mpcs2,azels,range,data.prm.cont,ti)
    case 'pmp',  PlotMpMapA(data.mpcc,data.mpcs,azels,range,data.prm.cont,ti)
    end
else
    switch data.prm.type
    case 'cmp1', PlotMpMapS(data.mpcc1,data.mpcs1,azels,range,data.prm.cont,ti)
    case 'cmp2', PlotMpMapS(data.mpcc2,data.mpcs2,azels,range,data.prm.cont,ti)
    case 'pmp',  PlotMpMapS(data.mpcc,data.mpcs,azels,range,data.prm.cont,ti)
    end
end

% plot multipath profile in skymap ---------------------------------------------
function PlotMpMapS(mpcc,mpcs,azels,range,cont,ti)
[fname,fsize]=gut('getfont','g');
clf, hold on
x=-90:3:90;
for n=1:length(x)
    for m=1:length(x)
        map(n,m)=rcvmpc([atan2(x(m),x(n)),(90-norm([x(n),x(m)]))*pi/180],mpcc,mpcs);
    end
end
step=[-100,range(1):(range(2)-range(1))/20:range(2),100];

[c,h]=contourf_v6(x,x,map,step); set(h,'edgecolor','none');
if cont
    [c,h]=contour_v6(x,x,map,step); set(h,'edgecolor',[0.5,0.5,0.5]);
    clabel(c,h,'fontname',fname,'fontsize',fsize); 
end
caxis(range), h=colorbar;
set(h,'position',[0.88,0.05,0.02,0.88],'fontname',fname,'fontsize',fsize)

set(gca,'position',[0.06,0.05,0.88,0.88],'fontname',fname,'fontsize',fsize)

for n=1:length(azels)
    if ~isempty(azels{n})
        i=find(abs(azels{n}(1:end-1,1)-azels{n}(2:end,1))>5*pi/180&...
               azels{n}(1:end-1,2)<30*pi/180);
        azels{n}(i,:)=nan;
        ggt('skyplot',azels{n},'-');
    end
end
ggt('skymap','noback','mask','fontname',fname,'fontsize',fsize);
title(ti)

% plot multipath profile as az/el graph ---------------------------------------
function PlotMpMapA(mpcc,mpcs,azels,range,cont,ti)
[fname,fsize]=gut('getfont','g');
clf, hold on
az=-180:5:180; el=0:5:90;
for n=1:length(el)
    for m=1:length(az)
        map(n,m)=rcvmpc([az(m)*pi/180,el(n)*pi/180],mpcc,mpcs);
    end
end
step=[-100,range(1):(range(2)-range(1))/20:range(2),100];

[c,h]=contourf_v6(az,el,map,step); set(h,'edgecolor','none');
if cont
    [c,h]=contour_v6(az,el,map,step); set(h,'edgecolor',[0.5,0.5,0.5]);
    clabel(c,h,'fontname',fname,'fontsize',fsize); 
end
for el=15:15:75
    for az=-180:30:180, plot([az,az],[-90,90],'k:'), end
    plot([-180,180],[el,el],'k:')
end
plot([-180,180,180,-180,-180],[0,0,90,90,0],'k')
set(gca,'position',[0.08,0.07,0.84,0.85],'xlim',[-180,180],'ylim',[0,90],...
    'xtick',-180:30:180,'ytick',0:15:90,'fontname',fname,'fontsize',fsize)
xlabel('Azimuth (deg)'), ylabel('Elevation (deg)')

caxis(range), h=colorbar('horiz');
set(h,'position',[0.08,0.08,0.84,0.025],'fontname',fname,'fontsize',fsize)

for n=1:length(azels)
    if ~isempty(azels{n})
        az=azels{n}(:,1); el=azels{n}(:,2);
        plot(az*180/pi,el*180/pi,'.','markersize',1)
    end
end
title(ti)

% read multipath profile -------------------------------------------------------
function ReadData(file)
data=get(gcf,'userdata');
if ~exist(file), return, end
td=[]; time=[]; rcv=''; mpcc1=[]; mpcc2=[]; mpcs1=[]; mpcs2=[]; mpcc=[]; mpcs=[];
azels={};
load(file)
data.td   =td;
data.time =time;
data.rcv  =rcv;
data.mpcs1=mpcs1;
data.mpcs2=mpcs2;
data.mpcc1=mpcc1;
data.mpcc2=mpcc2;
data.mpcc =mpcc;
data.mpcs =mpcs;
data.azels=azels;
if ~isempty(mpcc), data.prm.type='pmp'; else data.prm.type='cmp1'; end
gut('unchkmenu','plot');
gut('chkmenu','plot',data.prm.type);
set(gcf,'userdata',data,'name',['Multipath Profile : ',file]);

% time string ------------------------------------------------------------------
function str=tstr(ep), str=sprintf('%04d/%02d/%02d %02d:%02d',ep(1:5));

% matlab 6 compatible polts ----------------------------------------------------
function [c,h]=contour_v6(varargin)
v=version; v=str2num(v(1));
if v<=6, [c,h]=contour(varargin{:}); else [c,h]=contour('v6',varargin{:}); end

function [c,h]=contourf_v6(varargin)
v=version; v=str2num(v(1));
if v<=6, [c,h]=contourf(varargin{:}); else [c,h]=contourf('v6',varargin{:}); end
