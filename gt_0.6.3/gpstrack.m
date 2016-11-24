function h=gpstrack(varargin)
%-------------------------------------------------------------------------------
% [system] : GpsTools
% [module] : plot satellite ground tracks and receiver positions
% [func]   : plot satellite ground tracks and receiver positions
% [argin]  : 'prm',prm     : parameters
% [argout] : h = figure handle
% [note]   :
% [version]: $Revision: 32 $ $Date: 06/07/16 11:25 $
%            Copyright(c) 2004-2006 by T.Takasu, all rights reserved
% [history]: 04/11/18   0.1  new
%            06/03/06   0.2  add satellite visibility
%-------------------------------------------------------------------------------
if nargin>0&strncmp(varargin{1},'cb',2), feval(varargin{:})
else h=PlotMain(varargin{:}); end

% main -------------------------------------------------------------------------
function h=PlotMain(varargin)
h=gut('newfigg','','Satellite Track/Receiver Positions',[600,400,0,-68],...
      [mfilename,' cbClose']);
if isempty(h), return, end
prm=loadprm('prm_gpstrack','prm_gpstrack_def');
if nargin>=2&strcmp(varargin{1},'prm'), prm=varargin{2}; end
cb1=[mfilename,' cbPlot1'];
cb2=[mfilename,' cbPlot2'];
gut('newmenu','data','&Data ',...
    {'&Time Span...','S&atellites...','R&eceivers...','&New Window','-',...
     '&Export Plot...','-','&Map Area...','&Options...'},{},...
    {[mfilename,' cbTime'],[mfilename,' cbSelSat'],[mfilename,' cbSelRcv'],...
     [mfilename,' cbNew'],'',[mfilename,' cbExport'],'', [mfilename,' cbMap'],...
     [mfilename,' cbOpt']});
gut('newmenu','plot','&Plot',...
    {'Satellite &Ground Tracks','Receiver &Positions','Coverage &Areas',...
     '&Baseline Length','&Voronoi Region','-','Receiver &Sky Plot',...
     'Satellite &Visibility','&Elevation Angle by Stas','Ele&vation Angle by Sats',...
     'Satellite Eclipse &Period','Satellite-&Sun Angle','DOP'},...
    {'track','pos','area','tri','voro','','rcvsky','rcvvis','rcvel',...
     'satel','ecl','satsun','dop'},...
    {cb1,cb1,cb1,cb1,cb1,'',cb2,cb2,cb2,cb2,cb2,cb2,cb2});
gut('newmenu','sats','&Satellite',{'ALL',prm.sats{:}},{'ALL',prm.sats{:}},[mfilename,' cbSat']);
gut('newmenu','rcvs','&Receiver',{'ALL',prm.rcvs{:}},{'ALL',prm.rcvs{:}},[mfilename,' cbRcv']);
gut('chkmenu','plot',prm.type);
gut('chkmenu','sats',prm.sat);
gut('chkmenu','rcvs',prm.rcv);
set(gcf,'userdata',prm); UpdatePlot;

% callback on new --------------------------------------------------------------
function cbNew
feval(mfilename,'prm',get(gcf,'userdata'));

% callback on close ------------------------------------------------------------
function cbClose
SaveSetting;
closereq

% callback on menu time --------------------------------------------------------
function cbTime
prm=get(gcf,'userdata');
prm1={
't','Start Time (GPST)',   prm.tstart
't','End Time (GPST)',     prm.tend
'n','Time Interval (sec)', prm.tint
};
gut('newdlg','','Time Span',[325,120]);
gut('newprms','prm1',[15,35,300,27],prm1);
gut('newokcancelbtn','',[137,4,180,23]);
if ~gut('waitok'), return, end
[prm.tstart,prm.tend,prm.tint]=gut('getprms','prm1');
close
set(gcf,'userdata',prm); UpdatePlot;

% callback on select satellites ------------------------------------------------
function cbSelSat
prm=get(gcf,'userdata');
[prm.sats,ok]=editlist('Satellites',prm.sats);
if ~ok, return, end
delete(gut('geth','sats'));
gut('newmenu','sats','&Satellite',{'ALL',prm.sats{:}},{'ALL',prm.sats{:}},[mfilename,' cbSat']);
prm.sat='ALL';
gut('chkmenu','sats',prm.sat);
set(gcf,'userdata',prm);
UpdatePlot

% callback on select receiver --------------------------------------------------
function cbSelRcv
prm=get(gcf,'userdata');
[prm.rcvs,ok]=editlist('Receivers',prm.rcvs,'rcv');
if ~ok, return, end
delete(gut('geth','rcvs'));
gut('newmenu','rcvs','&Receivers',{'ALL',prm.rcvs{:}},{'ALL',prm.rcvs{:}},[mfilename,' cbRcv']);
prm.rcv='ALL';
gut('chkmenu','rcvs',prm.rcv);
set(gcf,'userdata',prm);
UpdatePlot

% callback on menu export plot -------------------------------------------------
function cbExport
persistent file dirs
if isempty(file), file='gpstrack'; end, if isempty(dirs), dirs=''; end
f=gcf; prm=get(f,'userdata');
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
if ~isempty(strmatch('sat',prm.type))&strcmp(prm.sat,'ALL')
    sat=prm.sat;
    for n=1:length(prm.sats)
        prm.sat=prm.sats{n}; data.prm=prm; set(f,'userdata',data)
        figure(f), UpdatePlot;
        print(format{:},opts{:},sprintf('%s%02d',fs,n));
    end
    prm.sat=sat; data.prm=prm; set(f,'userdata',data); UpdatePlot;
elseif ~isempty(strmatch('rcv',prm.type))&strcmp(prm.rcv,'ALL')
    rcv=prm.rcv;
    for n=1:length(prm.rcvs)
        prm.rcv=prm.rcvs{n}; set(f,'userdata',prm)
        figure(f), UpdatePlot;
        print(format{:},opts{:},sprintf('%s%02d',fs,n));
    end
    prm.rcv=rcv; set(f,'userdata',prm); UpdatePlot;
else
    print(format{:},opts{:},fs);
end

% callback on menu map area ----------------------------------------------------
function cbMap
prm=get(gcf,'userdata');
[prm.map,ok]=editmap('Map Area',prm.map);
if ok, set(gcf,'userdata',prm), UpdatePlot, end

% callback on menu options -----------------------------------------------------
function cbOpt
sel1={{'None','Start Position','End Position'},0:2};
sel2={{'OFF','Code Name','Full Name','Both Name'},0:3};
sel3={{'Dot','Line','Line and Dot'},1:3};
sel5={{'Combined','RINEX','RINEX(3H)','RINEX(1H)','RINEX(15min)'},...
      {'brdc','rinex','rinex3','rinex1','rinexh'}};
prm=get(gcf,'userdata');
prm1={
'n','Min Elevation Angle (deg)',     10
'p','Show Satellite Position in Map',sel1
'p','Show Receivers Name',           sel2
'p','Plot Line Type',                sel3
's','Plot Line Width / Marker Size', [0.5,3]
'c','Receiver Position Color',       'r'
'g','Receiver Name Font',            ''
'p','Input Navigation Message',      sel5
};
prm2={
' ','Navigation Message Directory',''
'd','',''
};
q='';
gut('newdlg','','Options',[320,262]);
gut('newprms','prm1',[12,72,300,23,120],prm1);
gut('newprms','prm2',[12,33,300,18,120],prm2);
gut('setprms','prm1',prm.elmin,prm.showf(1),prm.showf(2),prm.ptype,...
    prm.psize,prm.rcolor,prm.rfont,prm.src.nav);
gut('setprms','prm2',q,prm.dirs.nav,q);
gut('newokcancelbtn','',[133,4,180,23]);
if ~gut('waitok'), return, end
[prm.elmin,prm.showf(1),prm.showf(2),prm.ptype,prm.psize,prm.rcolor,prm.rfont,...
 prm.src.nav]=gut('getprms','prm1');
[q,prm.dirs.nav]=gut('getprms','prm2');
close
set(gcf,'userdata',prm);
SaveSetting; UpdatePlot;

% callback on menu plot --------------------------------------------------------
function cbPlot1
prm=get(gcf,'userdata');
gut('unchkmenu','plot',{'rcvsky','rcvvis','rcvel','satel','ecl','satsun','dop'});
gut('togglechk',gcbo);
prm.type=gut('getchk','plot');
set(gcf,'userdata',prm)
UpdatePlot;

function cbPlot2
prm=get(gcf,'userdata');
prm.type={get(gcbo,'userdata')};
gut('unchkmenu','plot');
gut('chkmenu','plot',prm.type);
set(gcf,'userdata',prm)
UpdatePlot;

% callback on menu satellite ---------------------------------------------------
function cbSat
prm=get(gcf,'userdata');
prm.sat=get(gcbo,'userdata');
gut('unchkmenu','sats');
gut('chkmenu','sats',prm.sat);
set(gcf,'userdata',prm);
UpdatePlot;

% callback on menu receiver ----------------------------------------------------
function cbRcv
prm=get(gcf,'userdata');
prm.rcv=get(gcbo,'userdata');
gut('unchkmenu','rcvs');
gut('chkmenu','rcvs',prm.rcv);
set(gcf,'userdata',prm);
UpdatePlot;

% save settingsrs --------------------------------------------------------------
function SaveSetting
prm=get(gcf,'userdata');
[path,file]=fileparts(which(mfilename));
save(fullfile(path,'settings','prm_gpstrack.mat'),'prm');

% update satellite ground track and receivers ----------------------------------
function UpdatePlot
prm=get(gcf,'userdata');
[td,ts]=caltomjd(prm.tstart); [tn,te]=caltomjd(prm.tend);
time=ts:prm.tint:te+(tn-td)*86400; if isempty(time), return, end
if strcmp(prm.sat,'ALL'), sats=prm.sats; else sats={prm.sat}; end
if strcmp(prm.rcv,'ALL'), rcvs=prm.rcvs; else rcvs={prm.rcv}; end
h=gcf; set(h,'pointer','watch')
if any(strcmp(prm.type,'rcvsky'))
    PlotSatSky(td,time,sats,rcvs,prm.elmin,prm.dirs,prm.src.nav,prm.showf,...
               prm.ptype,prm.psize);

elseif any(strcmp(prm.type,'rcvvis'))
    PlotStaVis(td,time,sats,rcvs,-32,prm.elmin,prm.dirs,prm.src.nav,prm.showf);

elseif any(strcmp(prm.type,'rcvel'))
    PlotStaEl(td,time,sats,rcvs,prm.elmin,prm.dirs,prm.src.nav,prm.showf);

elseif any(strcmp(prm.type,'satel'))
    PlotSatEl(td,time,sats,rcvs,prm.elmin,prm.dirs,prm.src.nav,prm.showf);

elseif any(strcmp(prm.type,'ecl'))
    PlotSatEcl(td,time,sats,rcvs,-32,prm.dirs,prm.src.nav);

elseif any(strcmp(prm.type,'satsun'))
    PlotSunAngle(td,time,sats,rcvs,-32,prm.dirs,prm.src.nav,prm.ptype);

elseif any(strcmp(prm.type,'dop'))
    if strcmp(prm.rcv,'ALL')&~isempty(prm.rcvs), rcv=prm.rcvs{1};
    else rcv=prm.rcv; end
    PlotDop(td,time,sats,rcv,-32,prm.dirs,prm.src.nav,prm.elmin);

else PlotTrack; end
set(h,'pointer','arrow')

% plot satellite ground track and receivers ------------------------------------
function PlotTrack
prm=get(gcf,'userdata');
[td,ts]=caltomjd(prm.tstart); [tn,te]=caltomjd(prm.tend);
time=ts:prm.tint:te+(tn-td)*86400; if isempty(time), return, end
if strcmp(prm.sat,'ALL'), sats=prm.sats; else sats={prm.sat}; end
if strcmp(prm.rcv,'ALL'), rcvs=prm.rcvs; else rcvs={prm.rcv}; end
clf, axis off
[fname,fsize]=gut('getfont','g');
p=get(gcf,'position'); pos=[0.001,0.001,0.998,(p(4)-20)/p(4)];
gmt('mmap','proj',prm.map.proj,'cent',prm.map.cent,'base',prm.map.base,...
    'scale',prm.map.scale,'pos',pos,'fontname',fname,'fontsize',fsize);
gmt('mcoast','lcolor',prm.map.color{1},'scolor',prm.map.color{2},'ccolor',...
    prm.map.color{3});
gmt('mgrid','gint',prm.map.gint(1),'lint',prm.map.gint(2),'color',prm.map.color{4});
ti=''; gpos=[];
posr=RcvPos(td,time(1),rcvs);
for n=1:length(rcvs), gpos(n,:)=eceftogeod(posr(n,:)'); end
if any(strcmp(prm.type,'area'))
    drawcover(rcvs,gpos,prm.elmin),
    for n=1:length(rcvs), drawarea(gpos(n,1),gpos(n,2),prm.elmin,prm.psize(1)); end
    ti=['Coverage Areas (Min Elev=',num2str(prm.elmin),'deg)'];
end
if any(strcmp(prm.type,'tri'))&~isempty(gpos)
    xl=get(gca,'xlim'); yl=get(gca,'ylim');
    [x,y,z]=gmt('lltoxy',gpos(:,2),gpos(:,1));
    i=find(~isnan(x)&~isnan(y));
    if length(i)==2
        plot(x,y,'-','color',prm.rcolor);
        d=norm(geodtoecef(gpos(1,:))-geodtoecef(gpos(2,:)));
        text(mean(x),mean(y),sprintf('%.1f',d/1E3),'fontname',fname,'fontsize',fsize,...
             'horizontal','center','vertical','middle');
    elseif length(i)>2
        for t=delaunay(x(i),y(i))'
            t=[t;t(1)]; plot(x(t),y(t),'-','color',prm.rcolor);
            for i=[1,2,3;2,3,1]
                if all(xl(1)<=x(t(i))&x(t(i))<=xl(2)&yl(1)<=y(t(i))&y(t(i))<=yl(2))
                    d=norm(geodtoecef(gpos(t(i(1)),:))-geodtoecef(gpos(t(i(2)),:)));
                    text(mean(x(t(i))),mean(y(t(i))),sprintf('%.1f',d/1E3),...
                         'fontname',fname,'fontsize',fsize,'horizontal','center',...
                         'vertical','middle');
                end
            end
        end
    end
end
if any(strcmp(prm.type,'voro'))&~isempty(gpos)
    xl=get(gca,'xlim'); yl=get(gca,'ylim');
    [x,y,z]=gmt('lltoxy',gpos(:,2),gpos(:,1));
    i=find(xl(1)<=x&x<=xl(2)&yl(1)<=y&y<=yl(2)&z>=0);
    if length(i)>=3
        h=voronoi(x(i)',y(i)','-'); set(h,'color',[0.7,0.7,0.7]);
    end
end
if any(strcmp(prm.type,'pos'))&~isempty(gpos)
    xl=get(gca,'xlim'); yl=get(gca,'ylim');
    [x,y,z]=gmt('lltoxy',gpos(:,2),gpos(:,1));
    i=find(xl(1)<=x&x<=xl(2)&yl(1)<=y&y<=yl(2)&z>=0);
    plot(x(i),y(i),'.','color',prm.rcolor,'markersize',prm.psize(2));
    for n=i'
        if any(prm.showf(2)==[1,3]), rcv=rcvs{n}; else rcv=''; end
        if any(prm.showf(2)==[2,3]), rcv=[rcv,' ',RcvName(rcvs{n})]; end
        if prm.showf(2)>0
            h=gmt('mtext',gpos(n,2),gpos(n,1),rcv,'horizontal','center',...
                'vertical','top');
            if isstruct(prm.rfont), set(h,prm.rfont); end
        end
    end
    if ~isempty(ti), ti=[' / ',ti]; end
    ti=['Receiver Positions',ti];
end
if any(strcmp(prm.type,'track'))
    poss=SatPos(td,time,sats,rcvs,prm.dirs,prm.src.nav);
    cs='bgrcmk'; h=[]; label={};
    for n=1:length(sats)
        lat=repmat(nan,length(time),1); lon=lat;
        for m=1:length(time)
            gpos=eceftogeod(poss(m,:,n)');
            lat(m)=gpos(1);
            lon(m)=gpos(2);
        end
        if all(~isnan(lon))
            c=cs(mod(n-1,6)+1);
            if prm.ptype==1|prm.ptype==3, hh=gmt('mplot',lon,lat,c,'marker','.','markersize',prm.psize(2)); end
            if prm.ptype==2|prm.ptype==3, hh=gmt('mplot',lon,lat,c,'linestyle','-','linewidth',prm.psize(1)); end
            if ~isempty(hh), h=[h;hh]; label={label{:},sats{n}}; end
            if prm.showf(1)>0
                if prm.showf(1)==1, lons=lon(1); lats=lat(1);
                else lons=lon(end); lats=lat(end); end
                gmt('mplot',lons,lats,c,'marker','.','markersize',prm.psize(2)*3);
                gmt('mtext',lons,lats,sats{n},'horizontal','center','vertical','top',...
                    'fontname',fname,'fontsize',fsize);
            end
        end
    end
    if ~isempty(h), h=legend(h,label); set(h,'color','w'); end
    ti=tstr(td,time(1)); if length(time)>1, ti=[ti,'-',tstr(td,time(end))]; end
    ti=[' Satellite Ground Tracks : ',ti];
end
title(ti);

% plot receiver skyplot --------------------------------------------------------
function PlotSatSky(td,time,sats,rcvs,elmin,dirs,nav,showf,ptype,psize)
[fname,fsize]=gut('getfont','g');
clf, set(gca,'position',[0.06,0.05,0.88,0.88],'fontname',fname,'fontsize',fsize)
ggt('skymap');
if isempty(rcvs), rcv=''; else rcv=rcvs{1}; end
posr=RcvPos(td,time(1),{rcv})';
poss=SatPos(td,time,sats,rcv,dirs,nav);
cs='bgrcmk'; h=[]; ns=0; label={};
for n=1:length(sats)
    if ~isnan(posr)
        azel=repmat(nan,length(time),2);
        for m=1:length(time), azel(m,:)=satazel(poss(m,:,n)',posr); end
        azel(azel(:,2)<0)=nan;
        if any(~isnan(azel))
            c=cs(mod(n-1,6)+1); ns=ns+1; label={label{:},sats{n}};
            if ptype==1|ptype==3, h(ns)=ggt('skyplot',azel,[c,'.'],'markersize',psize(2)); end
            if ptype==2|ptype==3, h(ns)=ggt('skyplot',azel,[c,'-'],'linewidth',psize(1)); end
            if showf(1)>0
                if showf(1)==1, ae=azel(1,:); else ae=azel(end,:); end
                ggt('skyplot',ae,[c,'.'],'markersize',psize(2)*3);
                ggt('skytext',ae,sats{n},'horizontal','center','vertical','top');
            end
        end
    end
end
if ~isempty(h)
    h=legend(h,label); pos=get(h,'position');
    set(h,'position',[[0.97,0.95]-pos(3:4),pos(3:4)]);
end
if showf(2)==2, rcv=RcvName(rcv); end
PlotTitle(['Sky Plot ',rcv],td,time)

% plot receiver visivility -----------------------------------------------------
function PlotStaVis(td,time,sats,rcvs,utc_tai,elmin,dirs,nav,showf)
clf
margin=[0.08,0.03,0.04,0.012];
ggt('subplotv','',1,1,'move',[1,0,1,0],'link',[1,0,1,0],'taxis',td,...
    'xlim',time([1,end])/3600,'margin',margin);
if isempty(rcvs), rcv=''; else rcv=rcvs{1}; end
posr=RcvPos(td,time(1),{rcv})';
poss=SatPos(td,time,sats,rcv,dirs,nav);
nx=length(sats); label='';
for n=1:length(time)
    U=ecsftoecef(td+(time(n)+19+utc_tai)/86400);
    [rsun,rmoon]=sunmoonpos(td+(time(n)+19+utc_tai)/86400);
    for m=1:nx
        ecl(n,m)=~isnan(poss(n,1,m))&shadowfunc(poss(n,:,m)',U*rsun,U*rmoon)<1;
    end
end
nobs=zeros(length(time),1);
for n=1:nx
    vis=repmat(nan,length(time),1);
    for m=1:length(time)
        azel=satazel(poss(m,:,n)',posr);
        if azel(2)*180/pi>=elmin, vis(m)=0; nobs(m)=nobs(m)+1; end
    end
    plot(time/3600,90*(nx-n)+45+vis,'b-','linewidth',3)
    vis(~ecl(:,n))=nan;
    plot(time/3600,90*(nx-n)+45+vis,'m-','linewidth',3)
    plot([-1E6,1E6],[90,90]*(nx-n),'k:')
    label=[sats{n},'|',label];
end
if showf(2)==2, rcv=RcvName(rcv); end
if nx>0
    set(gca,'xgrid','on','ylim',[0,90*nx],'ygrid','off','ytick',45:90:nx*90,...
        'yticklabel',label,'ticklength',[0,0]);
end
PlotTitle(sprintf('Satellite Visibilities %s (Min Elev=%.0fdeg)',rcv,...
          elmin),td,time);
ggt('subplotv','',1,1,'move',[1,0,1,0],'link',[1,0,1,0],'taxis',td,...
    'topts','nolabel','margin',margin,'xlim',time([1,end])/3600,'ylim',[0,16],...
    'box','off','color','none','yaxislocation','right','ticklength',[1E-3,0])
grid off;
plot(time/3600,nobs,'r-');

% plot receiver elevation angles -----------------------------------------------
function PlotStaEl(td,time,sats,rcvs,elmin,dirs,nav,showf)
clf
margin=[0.08,0.03,0.04,0.012];
ggt('subplotv','',1,1,'move',[1,0,1,0],'link',[1,0,1,0],'taxis',td,...
    'xlim',time([1,end])/3600,'margin',margin);
if isempty(rcvs), rcv=''; else rcv=rcvs{1}; end
posr=RcvPos(td,time(1),{rcv})';
poss=SatPos(td,time,sats,rcv,dirs,nav);
label=''; nx=length(sats);
for n=1:nx
    if ~isnan(posr)
        el=repmat(nan,length(time),1);
        for m=1:length(time)
            azel=satazel(poss(m,:,n)',posr);
            if azel(2)*180/pi>=0, el(m)=azel(2)*180/pi; end
        end
        plot(time/3600,el+90*(nx-n),':')
        el(el<elmin)=nan;
        plot(time/3600,el+90*(nx-n),'-','linewidth',2)
    end
    plot([-1E6,1E6],[90,90]*(nx-n),'k:')
    label=[sats{n},'|',label];
end
if nx>0
    set(gca,'xgrid','on','ylim',[0,90*nx],'ygrid','off','ytick',45:90:nx*90,...
        'yticklabel',label,'ticklength',[0,0]);
end
if showf(2)==2, rcv=RcvName(rcv); end
PlotTitle(['Elevation Angle ',rcv],td,time);

% plot satellite elevation angles ----------------------------------------------
function PlotSatEl(td,time,sats,rcvs,elmin,dirs,nav,showf)
clf
margin=[0.08,0.03,0.04,0.012];
ggt('subplotv','',1,1,'move',[1,0,1,0],'link',[1,0,1,0],'taxis',td,...
    'xlim',time([1,end])/3600,'margin',margin);
if isempty(sats), sat=''; else sat=sats{1}; end
posr=RcvPos(td,time(1),rcvs);
poss=SatPos(td,time,sat,rcvs,dirs,nav);
nx=length(rcvs); label='';
for n=1:nx
    if ~isnan(posr(n,:))
        el=repmat(nan,length(time),1);
        if ~isempty(poss)
            for m=1:length(time)
                azel=satazel(poss(m,:)',posr(n,:)');
                if azel(2)*180/pi>=0, el(m)=azel(2)*180/pi; end
            end
            plot(time/3600,el+90*(nx-n),':')
            el(el<elmin)=nan;
            plot(time/3600,el+90*(nx-n),'-','linewidth',2)
        end
    end
    plot([-1E6,1E6],[90,90]*(nx-n),'k:')
    if showf(2)==2, rcv=RcvName(rcvs{n}); else rcv=rcvs{n}; end
    label=[rcv,'|',label];
end
if nx>0
    set(gca,'xgrid','on','ylim',[0,90*nx],'ygrid','off','ytick',45:90:nx*90,...
        'yticklabel',label,'ticklength',[0,0]);
end
PlotTitle(['Elevation Angle ',sat],td,time);

% plot satellite eclipse periods -----------------------------------------------
function PlotSatEcl(td,time,sats,rcvs,utc_tai,dirs,nav)
clf
margin=[0.08,0.03,0.04,0.012];
ggt('subplotv','',1,1,'move',[1,0,1,0],'link',[1,0,1,0],'taxis',td,...
    'xlim',time([1,end])/3600,'margin',margin);
poss=SatPos(td,time,sats,rcvs,dirs,nav);
label=''; nx=length(sats);
for n=1:nx
    ecl=repmat(nan,length(time),1);
    for m=1:length(time)
        U=ecsftoecef(td+(time(m)+19+utc_tai)/86400);
        [rsun,rmoon]=sunmoonpos(td+(time(m)+19+utc_tai)/86400);
        if ~isnan(poss(m,1,n))&shadowfunc(poss(m,:,n)',U*rsun,U*rmoon)<1
            ecl(m)=0;
        end
    end
    plot(time/3600,100*(nx-n)+50+ecl,'b-','linewidth',3)
    plot([-1E6,1E6],[100,100]*(nx-n),'k:')
    label=[sats{n},'|',label];
end
if nx>0
    set(gca,'xgrid','on','ylim',[0,100*nx],'ygrid','off','ytick',50:100:nx*100,...
        'yticklabel',label,'ticklength',[0,0]);
end
PlotTitle('Eclipse Periods',td,time);

% plot satellite sun angle (alpha,beta) ----------------------------------------
function PlotSunAngle(td,time,sats,rcvs,utc_tai,dirs,src,ptype)
clf
if isempty(sats), sat=''; else sat=sats{1}; end
[nav,inav]=readnav(td,time,sat,rcvs,dirs.nav,src);
sunang=zeros(length(time),2);
for m=1:length(time)
    [poss,dts,vels]=navtostate(td,time(m),nav);
    [U,P,N,gmst,dx,dy,du]=ecsftoecef(td+(time(m)+19+utc_tai)/86400);
    rsat=U'*poss;
    vsat=U'*vels+du'*poss;
    rsun=sunmoonpos(td+(time(m)+19+utc_tai)/86400);
    psun=OrbitPos(rsat,vsat)*rsun;
    sunang(m,:)=[atan2(psun(2),psun(1)),asin(psun(3)/norm(psun))]*180/pi;
end
margin=[0.08,0.03,0.04,0.012];
ggt('subplotv','',1,2,'move',[1,0,1,0],'link',[1,0,1,0],'taxis',td,'topts','nolabel',...
    'xlim',time([1,end])/3600,'margin',margin);
if ptype==1|ptype==3, plot(time/3600,sunang(:,1),'.'), end
if ptype==2|ptype==3, plot(time/3600,sunang(:,1),'-'), end
ylim([-180,180]), ylabel('Alpha Angle (deg)')
PlotTitle(['Sun Angle ',sat],td,time)

ggt('subplotv','',2,2,'move',[1,0,1,0],'link',[1,0,1,0],'taxis',td,...
    'xlim',[time(1),time(end)]/3600,'margin',margin);
if ptype==1|ptype==3, plot(time/3600,sunang(:,2),'.'), end
if ptype==2|ptype==3, plot(time/3600,sunang(:,2),'-'), end
ylim([-90,90]), ylabel('Beta Angle (deg)')

% plot dop --------------------------------------------------------------------
function PlotDop(td,time,sats,rcv,utc_tai,dirs,nav,elmin)
posr=RcvPos(td,time(1),rcv)';
[gpos,E]=eceftogeod(posr);
poss=SatPos(td,time,sats,rcv,dirs,nav);
dop=repmat(nan,length(time),4); ecl=zeros(length(time),length(sats));
for n=1:length(time)
    U=ecsftoecef(td+(time(n)+19+utc_tai)/86400);
    [rsun,rmoon]=sunmoonpos(td+(time(n)+19+utc_tai)/86400);
    for m=1:length(sats), ecl(n,m)=shadowfunc(poss(n,:,m)',U*rsun,U*rmoon)<1; end
end
for n=1:length(time)
    A=[];
    for m=1:length(sats)
        azel=satazel(poss(n,:,m)',posr);
        if azel(2)>=elmin*pi/180&~ecl(n,m)
            e=E*(poss(n,:,m)'-posr);
            A=[A;-e'/norm(e),1];
        end
    end
    if size(A,1)>=4
        q=diag(inv(A'*A));
        dop(n,1:3)=sqrt([sum(q(1:3)),sum(q(1:2)),q(3)]);
    end
    dop(n,4)=size(A,1);
end
label={'PDOP','HDOP','VDOP','No of Sats'};
clf
margin=[0.08,0.03,0.04,0.012];
ggt('subplotv','',1,1,'move',[1,0,1,0],'link',[1,0,1,0],'taxis',td,...
    'xlim',time([1,end])/3600,'ylim',[0,15],'margin',margin);
h=plot(time/3600,dop); legend(h,label);
PlotTitle(sprintf('Dilution of Precision %s (Min Elev=%.0fdeg)',rcv,elmin),td,time)
ylabel('DOP / No of Sats')

% satellite position -----------------------------------------------------------
function poss=SatPos(td,time,sats,rcvs,dirs,src)
[nav,inav]=readnav(td,time,sats,rcvs,dirs.nav,src);
poss=repmat(nan,[length(time),3,length(sats)]);
for n=1:length(sats)
    for m=1:length(time)
        navi=nav(find(inav==n),:);
        if ~isempty(navi), poss(m,:,n)=navtostate(td,time(m),navi); end
    end
end

% read receiver positions ------------------------------------------------------
function posr=RcvPos(td,ts,rcvs)
posr=readpos(td,ts,rcvs,'','approx');
if isempty(posr), posr=[]; else posr=shiftdim(posr,1)'; end

% ecsf->orbital coord. ---------------------------------------------------------
function E=OrbitPos(rsat,vsat)
ex=rsat/norm(rsat);
ez=cross(rsat,vsat); ez=ez/norm(ez);
ey=cross(ez,ex);
E=[ex';ey';ez'];

% draw coverage boundaries -----------------------------------------------------
function drawarea(lat,lon,elmin,lw)
elmin=elmin*pi/180; re=6371; h=20200;
az=0.01:0.01:2*pi;
g=pi/2-elmin-asin(re/(re+h)*cos(elmin));
lat=lat*pi/180; lon=lon*pi/180;
latp=zeros(length(az),1); lonp=latp;
for n=1:length(az)
    latp(n)=asin(cos(g)*sin(lat)+sin(g)*cos(lat)*cos(az(n)));
    dlon=acos((cos(g)-sin(lat)*sin(latp(n)))/(cos(lat)*cos(latp(n))));
    if az(n)<=pi, lonp(n)=lon+dlon; else lonp(n)=lon-dlon; end
end
gmt('mplot',lonp*180/pi,latp*180/pi,[0.8,0.8,0.8],'linewidth',0.5);

% draw coverage areas ----------------------------------------------------------
function drawcover(rcvs,gpos,elmin)
hh=20200000;         % GPS satellite altitude(m)
elmin=elmin*pi/180;
cs=get(gcf,'colormap'); cs=cs(39:2:end,:);
[xl,yl]=gmt('getlim'); xi=(xl(2)-xl(1))/100;
for n=1:size(gpos,1), epos(:,n)=geodtoecef(gpos(n,:)); end
for x=xl(1)+xi/2:xi:xl(2)-xi/2
for y=yl(1)+xi/2:xi:yl(2)-xi/2
    [lon,lat]=gmt('xytoll',x,y);
    if ~isnan(lon)
        n=0;
        for m=1:length(rcvs)
            azel=satazel(geodtoecef([lat,lon,hh]),epos(:,m));
            if azel(2)>=elmin&n<size(cs,1), n=n+1; end
        end
        if n>0, plot(x,y,'color',cs(n,:),'marker','x','markersize',5); end
    end
end
end
for n=1:size(cs,1)
    h(n)=plot(0,0,'color',cs(n,:),'linestyle','none','marker','x',...
              'markersize',6,'visible','off');
    if n<size(cs,1), label{n}=sprintf('# of rcvs : %d',n);
    else label{n}=sprintf('# of rcvs >%d',n-1); end
end
legend(h,label)

% receiver name ----------------------------------------------------------------
function name=RcvName(rcv)
prmr=prm_gpsrcvs;
i=min(find(strcmp(rcv,prmr(:,1))));
if isempty(i), name=rcv; else name=prmr{i,3}; end

% plot title -------------------------------------------------------------------
function PlotTitle(ti,td,time)
title([ti,' : ',tstr(td,time(1)),'-',tstr(td,time(end))]);

% date/time string -------------------------------------------------------------
function str=tstr(td,ts)
str=sprintf('%04d/%02d/%02d %02d:%02d:%02.0f',mjdtocal(td,ts));

