function h=plotleo(varargin)
%-------------------------------------------------------------------------------
% [system] : GpsTools
% [module] : display/plot leo satellite orbit estimation results
% [func]   : display/plot leo satellite orbit estimation results
% [argin]  : ('prm',prm) = parameters
% [argout] : h           = figurea handle
% [note]   :
% [version]: $Revision: 8 $ $Date: 06/07/13 17:46 $
%            Copyright(c) 2004-2006 by T.Takasu, all rights reserved
% [history]: 05/11/04  0.1  new
%-------------------------------------------------------------------------------
if nargin>0&strncmp(varargin{1},'cb',2), feval(varargin{:});
else h=PlotMain(varargin{:}); end

% main -------------------------------------------------------------------------
function h=PlotMain(varargin)
h=gut('newfigg','','LEO Satellite Orbit',[600,400,0,-68],[mfilename,' cbClose']);
if isempty(h), return, end
prm=loadprm('prm_plotleo','prm_plotleo_def');
if nargin>=2&strcmp(varargin{1},'prm'), prm=varargin{2}; end
cb2=[mfilename,' cbPlot2'];
gut('newmenu','data','&Data ',{'&Read...','&New Window','-','Map Area...',...
    '&Options...'},{},...
    {[mfilename,' cbRead'],[mfilename,' cbNew'],'',[mfilename,' cbMap'],...
     [mfilename,' cbOpt']});
gut('newmenu','plot','&Plot',{'Position Error &3D/R/A/C','Position on Map'},...
    {'poserr','posmap'},[mfilename,' cbPlot']);
gut('chkmenu','plot',prm.type);
data.pos=[]; data.err=[]; data.cov=[];
data.prm=prm; set(h,'userdata',data);
if length(varargin)>0, ReadData; UpdatePlot; end

% callback on new --------------------------------------------------------------
function cbNew
data=get(gcf,'userdata');
feval(mfilename,'prm',data.prm);

% callback on close ------------------------------------------------------------
function cbClose
SaveSetting;
closereq;

% callback on menu read --------------------------------------------------------
function cbRead
sel1={{'Estimated(Forward)','Estimated(Backward)','Estimated(Smoothed)'},...
      {'posf','posb','posfb'}};
sel2={{'GRACE A','GRACE B','CHAMP'},{'GRAA','GRAB','CHAM'}};
sel3={{'OFF','Interp','Smoothed'},0:2};
data=get(gcf,'userdata'); prm=data.prm;
prm1={
't','Start Time (GPST)',[]
't','End Time (GPST)',  []
'n','Time Interval',    300
'n','Processing Unit Time (hr)',  24
'p','Satellite',        sel2
'p','Satellite Position',sel1
'p','Interpolation/Smoothing',sel3
};
prm2={
' ','Satellite Position Directory',''
'd','',prm.dirs.est
' ','Satellite Ephemeris Directory',''
'd','',prm.dirs.eph
};
gut('newdlg','','Read Data',[314,294]);
gut('newprms','prm1',[12,114,295,25,118],prm1);
gut('newprms','prm2',[12,34,295,18,118],prm2);
gut('newokcancelbtn','',[127,4,180,23]);
gut('setprms','prm1',prm.tstart,prm.tend,prm.tint,prm.tunit,prm.sat,prm.fb,prm.interp);
if ~gut('waitok'), return, end
[prm.tstart,prm.tend,prm.tint,prm.tunit,prm.sat,prm.fb,prm.interp]=gut('getprms','prm1');
[q,prm.dirs.est,q,prm.dirs.eph]=gut('getprms','prm2');
close
data.prm=prm; set(gcf,'userdata',data)
SaveSetting; ReadData; UpdatePlot;

% callback on menu options -----------------------------------------------------
function cbOpt
sel1={{'OFF','ON'},0:1};
sel2={{'Dot','Line','Line and Dot'},1:3};
data=get(gcf,'userdata'); prm=data.prm;
prm1={
'p','Show Mean and RMS Error',        sel1
'n','Show Est. Std Dev (sigma,0:OFF)', 0
'p','Plot Type',                      sel2
's','Plot Line Width / Marker Size',  [0.1,10]
's','Offset Radial (m)',               0
's','Offset Along-Trk (m)',            0
's','Offset Cross-Trk (m)',            0
's','Plot Axis Range (m)',            [-0.2,0.2]
};
gut('newdlg','','Options',[314,240]);
gut('newprms','prm1',[12,35,294,25,100],prm1);
gut('newokcancelbtn','',[126,4,180,23]);
gut('setprms','prm1',prm.showf(1),prm.nsig,prm.ptype,prm.psize,prm.off(1),...
    prm.off(2),prm.off(3),prm.range);
if ~gut('waitok'), return, end
[prm.showf(1),prm.nsig,prm.ptype,prm.psize,prm.off(1),prm.off(2),prm.off(3),...
 range]=gut('getprms','prm1');
if range(1)<range(2), prm.range=range; else prm.range=[-0.2,0.2]; end
close
data.prm=prm; set(gcf,'userdata',data);
SaveSetting; UpdatePlot;

% callback on menu map area ----------------------------------------------------
function cbMap
data=get(gcf,'userdata'); prm=data.prm;
[prm.map,ok]=editmap('Map Area',prm.map);
if ok, data.prm=prm; set(gcf,'userdata',data); UpdatePlot; end

% callback on menu plot --------------------------------------------------------
function cbPlot
data=get(gcf,'userdata'); prm=data.prm;
gut('unchkmenu','plot');
set(gcbo,'checked','on');
prm.type=get(gcbo,'userdata');
data.prm=prm; set(gcf,'userdata',data);
UpdatePlot;

% save settingsrs --------------------------------------------------------------
function SaveSetting
data=get(gcf,'userdata'); prm=data.prm;
[path,file]=fileparts(which(mfilename));
save(fullfile(path,'settings','prm_plotleo.mat'),'prm');

% update plot ------------------------------------------------------------------
function UpdatePlot
data=get(gcf,'userdata'); prm=data.prm;
if strcmp(prm.type,'posmap'), UpdateMap; return, end
clf
if isempty(data.err), return, end
margin=[0.08,0.03,0.04,0.012];
ggt('subplotv','',1,2,'move',[1,1,1,1],'link',[1,0,1,0],'taxis',prm.td,...
    'topts','nolabel','xlim',[prm.time(1),prm.time(end)]/3600,...
    'ylim',[0,prm.range(2)*2],'margin',margin);
err3d=sqrt(sum(data.err.^2,2));
cov3d=sqrt(sum(data.cov.^2,2));
if prm.ptype==1|prm.ptype==3, plot(prm.time/3600,err3d,'.','markersize',prm.psize(2)); end
if prm.ptype==2|prm.ptype==3, plot(prm.time/3600,err3d,'-','linewidth',prm.psize(1)); end
if prm.nsig>0, plot(prm.time/3600,prm.nsig*sqrt(cov3d),'r:'); end
if prm.showf(1)
    ggt('mtext','_1',sprintf('REF: Level1B MEAN: %.4fm RMS: %.4fm',meann(err3d),rmsn(err3d)));
end
ylabel('Position Error 3D (m)');
ti=[tstr(prm.td,prm.time(1)),'-',tstr(prm.td,prm.time(end))];
title(['Satellite Orbit ',prm.sat,' : ', ti]);

ggt('subplotv','',2,2,'move',[1,1,1,1],'link',[1,0,1,0],'taxis',prm.td,...
    'xlim',[prm.time(1),prm.time(end)]/3600,'ylim',prm.range,'margin',margin);
if prm.ptype==1|prm.ptype==3, h=plot(prm.time/3600,data.err,'.','markersize',prm.psize(2)); end
if prm.ptype==2|prm.ptype==3, h=plot(prm.time/3600,data.err,'-','linewidth',prm.psize(1)); end
legend(h,{'Radial','Along-Track','Cross-Track'});
if prm.nsig>0
    plot(prm.time/3600,-prm.nsig*sqrt(data.cov),':');
    plot(prm.time/3600, prm.nsig*sqrt(data.cov),':');
end
if prm.showf(1)
    ggt('mtext','_2',sprintf('REF: Level1B\nMEAN R: %.4fm A: %.4fm C: %.4fm\nRMS R:%.4fm A:%.4fm C: %.4fm',...
        meann(data.err),rmsn(data.err)),3);
end
ylabel('Position Error R/A/C (m)');

% plot satellite position error map --------------------------------------------
function UpdateMap
data=get(gcf','userdata'); prm=data.prm;
clf
if isempty(data.pos), return, end
[fname,fsize]=gut('getfont','g');
p=get(gcf,'position'); pos=[0.001,0.001,0.998,(p(4)-20)/p(4)];
axis off, box on
gmt('mmap','proj',prm.map.proj,'cent',prm.map.cent,'base',prm.map.base,...
    'scale',prm.map.scale,'pos',pos,'fontname',fname,'fontsize',fsize);
gmt('mcoast','lcolor',prm.map.color{1},'scolor',prm.map.color{2},'ccolor',prm.map.color{3});
gmt('mgrid','gint',prm.map.gint(1),'lint',prm.map.gint(2),'color',prm.map.color{4});
[xl,yl]=gmt('getlim');
gpos=repmat(nan,size(data.pos));
for n=1:size(data.pos,1)
    gpos(n,:)=eceftogeod(data.pos(n,:)');
end
if prm.ptype==1|prm.ptype==3, gmt('mplot',gpos(:,2),gpos(:,1),prm.pcolor,'linestyle','none','marker','.','markersize',prm.psize(2)); end
if prm.ptype==2|prm.ptype==3, gmt('mplot',gpos(:,2),gpos(:,1),prm.pcolor,'linewidth',prm.psize(1)); end
ti=[tstr(prm.td,prm.time(1)),'-',tstr(prm.td,prm.time(end))];
title(['Satellite Orbit ',prm.sat,' : ',ti]);

% read data -------------------------------------------------------------------
function ReadData
data=get(gcf,'userdata'); prm=data.prm;
[prm.td,ts]=caltomjd(prm.tstart); [tn,te]=caltomjd(prm.tend);
prm.time=ts:prm.tint:te+(tn-prm.td)*86400; if isempty(prm.time), return, end
data.pos=[]; data.err=repmat(nan,length(prm.time),3); data.cov=data.err;
h=gcf; set(h,'pointer','watch');
[poss,covs]=readpos(prm.td,prm.time,prm.sat,prm.dirs.est,prm.fb,prm.tunit);
switch prm.interp
case 1
    i=find(~isnan(poss(:,1))&covs(:,1)<1);
    poss=interplag(prm.time(i),poss(i,:),prm.time,12);
otherwise
    i=find(covs(:,1)>=1); poss(i,:)=nan; covs(i,:)=nan;
end
t=[]; p=[];
switch prm.sat
case {'GRAA','GRAB'}
    tu=86400;
    for ts=(floor(prm.time(1)/tu):floor(prm.time(end)/tu))*tu
        ep=mjdtocal(prm.td+ts/86400);
        file=sprintf('GNV1B_%04d-%02d-%02d_%s_00.dat',ep(1:3),prm.sat(end));
        file=gfilepath(prm.dirs.eph,file,ep,prm.sat);
        [e,tt,k,pp,sig,s]=readgrace1b(file);
        if ~isempty(e), t=[t;tt+(caltomjd(e)-prm.td)*86400]; p=[p;pp]; end
    end
case 'CHAM'
    tu=43200; to=36000;
    for ts=(floor((prm.time(1)-to)/tu):floor((prm.time(end)-to)/tu))*tu+to
        ep=mjdtocal(prm.td,ts);
        doy=caltomjd(ep(1:3))-caltomjd([ep(1),1,1])+1;
        file=sprintf('CH-OG-3-RSO+CTS-CHA_%04d_%03d_%02d.dat',ep(1),doy,ep(4));
        file=gfilepath(prm.dirs.eph,file,ep,prm.sat);
        [e,tt,pp]=readchorb(file);
        if ~isempty(e)
            i=find(ep(4)*3600<=tt&tt<ep(4)*3600+tu);
            t=[t;tt(i)+(caltomjd(e)-prm.td)*86400]; p=[p;pp(i,:)];
        end
    end
end
set(h,'name',['LEO Satellite Orbit : ',prm.dirs.est,' (',prm.fb,')']);
set(h,'pointer','arrow'); figure(h);
if ~isempty(t)
    [tt,i,j]=intersect(prm.time,t);
    data.pos=poss; data.err=repmat(nan,length(prm.time),3);
    data.cov=data.err;
    err=poss(i,1:3)-p(j,1:3);
    for n=1:length(j)
        E=eceftosatf(p(j(n),:)');
        data.err(i(n),:)=err(n,:)*E'-prm.off;
        data.cov(i(n),:)=diag(E*diag(covs(n,:))*E')';
    end
end
data.prm=prm; set(gcf,'userdata',data);

% time str ---------------------------------------------------------------------
function s=tstr(td,ts)
ep=mjdtocal(td,ts); s=sprintf('%04d/%02d/%02d %02d:%02d',ep(1:5));

function r=meann(x), r=mean(x(~isnan(x(:,1)),:),1);
function r=rmsn(x), r=sqrt(mean(x(~isnan(x(:,1)),:).^2,1));
