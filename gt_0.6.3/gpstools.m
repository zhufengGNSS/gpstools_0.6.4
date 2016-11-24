function h=gpstools(varargin)
%-------------------------------------------------------------------------------
% [system] : GpsTools
% [module] : GpsTools main menu
% [func]   : GpsTools main menu
% [argin]  : option
%                'version' = get version
% [argout] : h = figure handle/version
% [note]   :
% [version]: $Revision: 27 $ $Date: 06/07/22 22:56 $
%            Copyright(c) 2004-2006 by T.Takasu, all rights reserved
% [history]: 04/11/15  0.1  new
%            06/03/06  0.2  supports GT0.6.1
%            06/07/08  0.3  supports GT0.6.3
%-------------------------------------------------------------------------------
ver='0.6.3'; % GpsTools version

if nargin<1, h=GpsToolsMain(ver);
elseif strcmp(varargin{1},'version') h=ver;
else feval(varargin{:}); end

% main menu --------------------------------------------------------------------
function h=GpsToolsMain(ver)
prm=loadprm('prm_gpstools','prm_gpstools_def');
h=gut('newfig','gtmain',['GpsTools ver.',ver],[600,23],[mfilename,' cbExit']);
if isempty(h), return, end
setgui(prm);
gut('newbtnh','menu',[1,1,600,22],...
    {'File...','Estimate...','Plot','Tools','Help','Exit'},...
    {[mfilename,' cbFile'],[mfilename,' cbEst'],[mfilename,' cbPlot'],...
     [mfilename,' cbTools'],[mfilename,' cbMisc'],[mfilename,' cbExit']})
gut('newcmenu','plot',[200,0],...
    {'&Observation Data...','-','&Satellite Orbit...','Receiver &Position...',...
     'Satellite/Receiver &Clock...','&Tropospheric Parameters...',...
     '&Grid Point Values...','&Ionospheric Parameters...',...
     '&LEO Satellite Orbit...','&Residuals/Statistics...','&Multipath Profile...',...
     '-','Satellite &Track/Receiver Pos...'},{},...
    {[mfilename,' cbPlotObs'],'',[mfilename,' cbPlotEph' ],[mfilename,' cbPlotPos'],...
     [mfilename,' cbPlotClk'],[mfilename,' cbPlotTrop'],[mfilename,' cbPlotGpv'],...
     [mfilename,' cbPlotIon'],[mfilename,' cbPlotLeo'],[mfilename,' cbPlotStats'],...
     [mfilename,' cbPlotMp'],'',[mfilename,' cbPlotTrack']});
gut('newcmenu','tools',[300,0],...
    {'&Generate Products...','Generate &PWV...','&Ocean Loading Parameters...',...
     'Generate &Multipath Profile...','-','&Download...'},{},...
    {[mfilename,' cbProd' ],[mfilename,' cbPwv'],[mfilename,' cbOload'],...
     [mfilename,' cbMultp'],'',[mfilename,' cbDown']});
gut('newcmenu','setting',[400,0],...
    {'&Settings...','&Help...','-','&About...'},{},...
    {[mfilename,' cbSetting'],[mfilename,' cbHelp'],'',[mfilename,' cbAbout']});
set(h,'userdata',prm);

% callback on main menu --------------------------------------------------------
function cbFile,      gpsfiles;
function cbEst,       execmd('ggpsest_');
function cbPlot,      gut('showmenu','plot');
function cbTools,     gut('showmenu','tools');
function cbMisc,      gut('showmenu','setting');

% callback on exit -------------------------------------------------------------
function cbExit
prm=get(gcf,'userdata'); prm.pos=get(gcf,'position');
[path,file]=fileparts(which(mfilename));
save(fullfile(path,'settings','prm_gpstools.mat'),'prm');
closereq

% callback on resize -----------------------------------------------------------
function cbResize
p=get(gcf,'position');
for n=1:6
    h=gut('geth',['menu_',num2str(n)]);
    set(h,'position',[p(3)/6*(n-1)+2,1,max(p(3)/6-2,50),max(p(4)-2,15)]);
end

% plot menu --------------------------------------------------------------------
function cbPlotObs,   plotobs;
function cbPlotEph,   ploteph;
function cbPlotPos,   plotpos;
function cbPlotClk,   plotclk;
function cbPlotTrop,  plottrop;
function cbPlotGpv,   plotgpv;
function cbPlotIon,   plotion;
function cbPlotLeo,   plotleo;
function cbPlotPwv,   plotpmap;
function cbPlotStats, plotstats;
function cbPlotMp,    plotmp;
function cbPlotTrack, gpstrack;

% generate product menu --------------------------------------------------------
function cbProd
prm=get(gcf,'userdata');
h=gut('newfig','gtprod','Generate Products',[314,308,0,-50],[mfilename,' cbProdClose']);
if isempty(h), return, end
set(gcf,'userdata',prm);
sel1={{'Sat Ephemeris (SP3)','Sat Clock(RINEX CLK)','Rcv Clock (RINEX CLK)',...
       'Rcv Position (SINEX)'},{'eph','clk','rclk','pos'}};
sel2={{'Est(Forward)','Est(Backward)','Est(Smoothed)'},{'f','b','fb'}};
sel3={{'','Include Velocity (SP3)'},{'','vel'}};
prm1={
't','Start Time (GPST)',           []
't','End Time (GPST)',             []
'n','Time Interval (sec)',         300
'n','Processing Unit Time (hr)',   24
'p','Product Type',                sel1
'p','Product Source',              sel2
'p','Generation Options',          sel3
'b','Satellites/Receivers',[mfilename,' cbSelSat']
};
prm2={
' ','Input Estimation Data Directory',''
'd','',prm.genprod.idir;
' ','Output Products Directory',''
'd','',prm.genprod.odir;
};
gut('newprms','prm1',[12,110,295,24,118],prm1);
gut('newprms','prm2',[12,32,295,18],prm2);
gut('newbtnh','btns', [127,4,180,22],{'Generate','Close'},{[mfilename,' cbProdExe'],''});
gut('setprms','prm1',prm.genprod.tstart,prm.genprod.tend,prm.genprod.tint,...
    prm.genprod.tunit,prm.genprod.type,prm.genprod.fb,prm.genprod.opt,'');
gut('setudata','prm1_8',prm.genprod.list);

function cbProdExe
prm=GetProdPrm;
[td,ts]=caltomjd(prm.tstart);
[tn,te]=caltomjd(prm.tend);
tspan=te-ts+(tn-td)*86400;
dspan=ceil(tspan/86400);

h=gcf; set(h,'pointer','watch'); gut('setdisall',h);
switch prm.type
case 'eph'
    ephtosp3('td',td,'span',dspan,'file','ephg','tint',prm.tint,...
             'idir',prm.idir,'odir',prm.odir,'fb',prm.fb,prm.opt);
case 'clk'
    dirs.est=prm.idir; dirs.clk=prm.odir;
    clktornx(prm.tstart,prm.tend,prm.tint,prm.tunit,dirs,['clk',prm.fb]);
case 'pos'
    postosinex(td,dspan,prm.tint,prm.list{2},prm.idir,prm.odir,...
               prm.fb,prm.tunit);
end
gut('setenaall',h); set(h,'pointer','arrow');

function cbProdClose
p=GetProdPrm;
closereq;
h=gut('geth','gtmain'); if isempty(h), return, end
prm=get(h,'userdata'); prm.genprod=p; set(h,'userdata',prm);

function prm=GetProdPrm
[prm.tstart,prm.tend,prm.tint,prm.tunit,prm.type,prm.fb,prm.opt,q]=gut('getprms','prm1');
[q,prm.idir,q,prm.odir]=gut('getprms','prm2');
prm.list=gut('getudata','prm1_8');

% callback on menu select satellite --------------------------------------------
function cbSelSat
list=get(gcbo,'userdata');
switch gut('getsel','prm1_5')
case {'eph','clk'},  [list{1},ok]=editlist('Satellites',list{1});
case {'rclk','pos'}, [list{2},ok]=editlist('Receivers',list{2},'rcv');
end
if ok, set(gcbo,'userdata',list); end

% generate ephemeris menu ------------------------------------------------------
function cbEph
prm=get(gcf,'userdata');
h=gut('newfig','gteph','Generate Ephemeris',[314,270,0,-50],[mfilename,' cbEphClose']);
if isempty(h), return, end
s1={'IGS Final','IGS Rapid','IGS URapid','IGS URapid(Pred)','CODE','EMR','ESA','GFZ','JPL','MIT','EPHG'};
c1={'igs','igr','igu','igp','cod','emr','esa','gfz','jpl','mit','ephg'};
s2={'','IERS BulltenB','IGS Final','IGS Rapid','IGS URapid'};
c2={'','bulb','igs','igr','igu'};
prm1={
'y','Ephemeris Start Day (GPST)', []
'y','Ephemeris End Day (GPST)',   []
'n','Ephemeris Interval (sec)',   30
'n','Ephemeris Unit Time (hr)',   24
'p','Input Ephemeris',       {s1,c1}
'p','Input ERP',             {s2,c2}
};
prm2={
' ','Input Ephemeris Directory',''
'd','',prm.geneph.ephdir
' ','Input ERP Directory',''
'd','',prm.geneph.erpdir
' ','Output Ephemeris Directory',''
'd','',prm.geneph.odir
};
gut('newprms','prm1',[12,110,295,24,118],prm1);
gut('newprms','prm2',[12,32,295,18],prm2);
gut('newbtnh','btns', [127,4,180,22],{'Generate','Close'},{[mfilename,' cbEphExe'],''});
gut('setprms','prm1',prm.geneph.tstart,prm.geneph.tend,prm.geneph.tint,...
    prm.geneph.tunit,prm.geneph.eph,prm.geneph.erp);

function cbEphExe
prm=GetEphPrm;
[td,ts]=caltomjd(prm.tstart); [tn,te]=caltomjd(prm.tend);
dirs.eph=prm.ephdir; dirs.erp=prm.erpdir;
h=gcf; set(h,'pointer','watch'); gut('setdisall',h);
geneph(td,'span',tn-td+1,'ephsrc',prm.eph,'erpsrc',prm.erp,...
       'dirs',dirs,'outdir',prm.odir,'tint',prm.tint,'tunit',prm.tunit);
gut('setenaall',h); set(h,'pointer','arrow');

function cbEphClose
p=GetEphPrm;
closereq
h=gut('geth','gtmain'); if isempty(h), return, end
prm=get(h,'userdata'); prm.geneph=p; set(h,'userdata',prm);

function prm=GetEphPrm
[prm.tstart,prm.tend,prm.tint,prm.tunit,prm.eph,prm.erp]=gut('getprms','prm1');
[q,prm.ephdir,q,prm.erpdir,q,prm.odir]=gut('getprms','prm2');

% generate pwv menu -----------------------------------------------------------
function cbPwv
prm=get(gcf,'userdata');
h=gut('newfig','gtpwv','Generate PWV',[314,266,0,-50],[mfilename,' cbPwvClose']);
if isempty(h), return, end
prm1={
't','Start Time (GPST)',  []
't','End Time (GPST)',    []
'n','Time Interval (sec)', 3600
'n','Processing Unit Time', 24
'b','Receivers', [mfilename,' cbSelRcv']
};
prm2={
' ','Input Estimation Data Directory',''
'd','',prm.genpwv.estdir
' ','Input Meteorological Data Directory',''
'd','',prm.genpwv.metdir
' ','Output Products Directory',''
'd','',prm.genpwv.odir
};
gut('newprms','prm1',[12,144,295,23,118],prm1);
gut('newprms','prm2',[12,32,295,18],prm2);
gut('newbtnh','btns', [127,4,180,22],{'Generate','Close'},...
    {[mfilename,' cbPwvExe'],[mfilename,' cbPwvClose']});
gut('setprms','prm1',prm.genpwv.tstart,prm.genpwv.tend,prm.genpwv.tint,...
    prm.genpwv.tunit,'');
gut('setudata','prm1_5',prm.genpwv.rcvs);

function cbPwvExe
prm=GetPwvPrm;
td=caltomjd(prm.time);
h=gcf; set(h,'pointer','watch'); gut('setdisall',h);
pwvgps(td,0:prm.tint:prm.span*86400-prm.tint,prm.rcvs,prm.dirs,prm.outdir,prm.tunit);
gut('setenaall',h); set(h,'pointer','arrow');

function cbPwvClose
p=GetPwvPrm;
closereq
h=gut('geth','gtmain'); if isempty(h), return, end
prm=get(h,'userdata'); prm.genpwv=p; set(h,'userdata',prm);

function cbSelRcv
rcvs=get(gcbo,'userdata');
[rcvs,ok]=editlist('Stations',rcvs,'rcv');
if ok, set(gcbo,'userdata',rcvs); end

function prm=GetPwvPrm
[prm.tstart,prm.tend,prm.tint,prm.tunit,q]=gut('getprms','prm1');
[q,prm.estdir,q,prm.metdir,q,prm.odir]=gut('getprms','prm2');
prm.rcvs=gut('getudata','prm1_5');

% ocean loading menu -----------------------------------------------------------
function cbOload
prm=get(gcf,'userdata');
h=gut('newfig','gtoload','Generate Ocean Loading Params',[300,340,0,-50],...
      [mfilename,' cbOloadClose']);
if isempty(h), return, end
sel1={{'Approx. Pos','ITRF2000','ITRF97','IGS Final','GSI Pos','Est(forward)',...
       'Est(backward)','Est(smoothed)'},...
      {'approx','itrf2000','itrf97','igssnx','gsipos','posf','posb','posfb'}};
prm1={
'p','Receiver Position', sel1
'y','Day of Position',   [2000,1,1];
'n','Processing Unit Time (hr)', 24
};
prm2={
' ','GOTIC2 Directory',''
'd','',prm.genoload.gdir
' ','Receiver Position Directory',''
'd','',prm.genoload.pdir
' ','Output Ocean Loading Parameter File',''
'f','',prm.genoload.file
};
gut('newtext','',[12,304,160,23],'Receivers');
gut('newbtn','',[92,310,78,20],'...',[mfilename,' cbOloadRcv']);
gut('newlist','rcvs',[177,213,119,120],{},{},2,'');
gut('newprms','prm1',[12,144,282,24,116],prm1);
gut('newprms','prm2',[12,32,282,18],prm2);
gut('newbtnh','btns', [114,4,180,22],{'Generate','Close'},{[mfilename,' cbOloadExe'],''});
gut('setprms','prm1',prm.genoload.psrc,prm.genoload.time,prm.genoload.tunit);
gut('setlist','rcvs',prm.genoload.rcvs);

function cbOloadRcv
rcvs=gut('getstring','rcvs');
gut('setlist','rcvs',editlist('Receivers',rcvs,'rcv'),{});

function cbOloadExe
prm=GetOloadPrm;
if isempty(prm.file)|isempty(prm.rcvs), return, end
[td,ts]=caltomjd(prm.time,0);
h=gcf; set(h,'pointer','watch'); gut('setdisall',h);
genblq(prm.file,td,ts,prm.rcvs,prm.pdir,prm.psrc,gfilepath('',prm.gdir),prm.tunit);
gut('setenaall',h); set(h,'pointer','arrow');

function cbOloadClose
p=GetOloadPrm;
closereq;
h=gut('geth','gtmain'); if isempty(h), return, end
prm=get(h,'userdata'); prm.genoload=p; set(h,'userdata',prm);

function prm=GetOloadPrm
[prm.psrc,prm.time,prm.tunit]=gut('getprms','prm1');
[q,prm.gdir,q,prm.pdir,q,prm.file]=gut('getprms','prm2');
prm.rcvs=gut('getstring','rcvs');

% multipath menu ---------------------------------------------------------------
function cbMultp
prm=get(gcf,'userdata');
h=gut('newfig','gtoload','Generate Multipath Profile',[314,322,0,-50],...
      [mfilename,' cbMultpClose']);
if isempty(h), return, end
sel1={{'Forward','Backward'},{'f','b'}};
prm1={
't','Start Time (GPST)',        []
't','End Time (GPST)',          []
'n','Time Interval (sec)',      300
'n','Processing Unit Time (H)', 24
'p','Estimation Direction',     sel1
'b','Stations', [mfilename,' cbMultpRcv']
'n','Degree of Spherical Harmonic', 12
};
prm2={
' ','Estimation Data Directory',''
'd','',prm.multp.estdir
' ','Navigation Message Directory',''
'd','',prm.multp.navdir
' ','Output Multipath Profile Directory',''
'd','',prm.multp.odir
};
gut('newprms','prm1',[12,148,295,24,118],prm1);
gut('newprms','prm2',[12,32,295,18],prm2);
gut('newbtnh','btns', [127,4,180,22],{'Execute','Close'},{[mfilename,' cbMultpExe'],''});
gut('setprms','prm1',prm.multp.tstart,prm.multp.tend,prm.multp.tint,prm.multp.tunit,...
    prm.multp.fb,'',prm.multp.nmax);
gut('setudata','prm1_6',prm.multp.rcvs);

function cbMultpRcv
rcvs=get(gcbo,'userdata');
rcvs=editlist('Stations',rcvs,'rcv');
set(gcbo,'userdata',rcvs);

function cbMultpExe
prm=GetMultpPrm;
[td,ts]=caltomjd(prm.tstart); [tn,te]=caltomjd(prm.tend);
time=ts:prm.tint:te+(tn-td)*86400-prm.tint;
dirs.est=prm.estdir; dirs.nav=prm.navdir;
h=gcf; set(h,'pointer','watch'); gut('setdisall',h);
gpsmpp(td,time,prm.rcvs,dirs,prm.odir,prm.fb,prm.tunit,prm.nmax);
gut('setenaall',h); set(h,'pointer','arrow');

function cbMultpClose
p=GetMultpPrm;
closereq;
h=gut('geth','gtmain'); if isempty(h), return, end
prm=get(h,'userdata'); prm.multp=p; set(h,'userdata',prm);

function prm=GetMultpPrm
[prm.tstart,prm.tend,prm.tint,prm.tunit,prm.fb,q,prm.nmax]=gut('getprms','prm1');
[q,prm.estdir,q,prm.navdir,q,prm.odir]=gut('getprms','prm2');
prm.rcvs=gut('getudata','prm1_6');

% download menu ----------------------------------------------------------------
function cbDown, execmd('gpsdown_');

% setting menu -----------------------------------------------------------------
function cbSetting
prm=get(gcf,'userdata');
prm1={
'g','UI Font',          {}
'g','Button Font',      {}
'g','Fixed Width Font', {}
'g','Graph Font',       {}
'n','Text Viewer Max Line Count', 20000
};
prm2={' ','Text Editor','';'f','','';' ','Help Browser','';'f','',''};
h=gut('newdlg','','Settings',[290,236,0,-148],1);
gut('newprms','prm1',[10,112,274,23],prm1);
gut('newprms','prm2',[10,32,274,18],prm2);
gut('newokcancelbtn','',[104,4,180,22]);
for n=1:4
    font{n}.FontName=prm.fonts{n}{1};
    font{n}.FontSize=prm.fonts{n}{2};
    font{n}.FontUnits='points';
    font{n}.FontWeight='normal';
    font{n}.FontAngle='normal';
end
gut('setprms','prm1',font{1},font{2},font{3},font{4},prm.maxtext);
gut('setprms','prm2','',prm.editor,'',prm.browser);
if ~gut('waitok'), return, end
[font{1},font{2},font{3},font{4},prm.maxtext]=gut('getprms','prm1');
[q,prm.editor,q,prm.browser]=gut('getprms','prm2');
for n=1:4, prm.fonts{n}={font{n}.FontName,font{n}.FontSize}; end
close
set(gcf,'userdata',prm);
setgui(prm);

function setgui(prm)
f={'u','b','f','g'}; 
for n=1:4, gut('setfont',f{n},prm.fonts{n}{1},prm.fonts{n}{2}); end
gutprm=gut('getgutprm');
gutprm.maxtext=prm.maxtext;
gutprm.editor=prm.editor;
gut('initgut',gutprm)

% help menu --------------------------------------------------------------------
function cbHelp
prm=get(gcf,'userdata');
[path,f,e]=fileparts(which(mfilename));
[s,w]=dos(['"',prm.browser,'" "',fullfile(path,'help','gpstools.htm'),'" &']);

% about menu -------------------------------------------------------------------
function cbAbout
[dirs,f]=fileparts(which(mfilename));
msg={get(gcf,'name'),'Copyright(c) 2004-2006, All rights reserved',...
     'T.Takasu and Kasai Design Office Ltd.',...
     '$Revision: 27 $ $Date: 06/07/22 22:56 $',...
     ['Path: ',dirs],['Matlab: ',version]};
gut('newdlg','','About...',[260,240],0);
gut('newokbtn','',[190,8,55,20]);
[fn,fs]=gut('getfont','u'); axis off
gmt('mmap','proj','ortho','cent',[135,40],'pos',[0.05,0.05,0.9,0.9]);
gmt('mcoast','ccolor',[0.7,0.7,0.7],'lcolor','none');
gmt('mgrid','gint',30,'lint',0,'color',[0.7,0.7,0.7]);
tr=SatTrack; gmt('mplot',tr(:,2),tr(:,1),[0,0.7,1])
[x,y]=gmt('lltoxy',135,40);
text(x,y,msg{1},'horizontal','center','fontname',fn,'fontweight','bold',...
    'fontsize',14);
for n=2:length(msg)
    text(x,y-0.125*n,msg{n},'fontname',fn,'fontsize',fs,'horizontal','center',...
         'interpreter','none');
end
if gut('waitok'), close, end

function tr=SatTrack
prm=loadprm('prm_gpsest'); [td,t]=caltomjd(prm.tstart); tr=[];
nav=readnav(td,0,sprintf('GPS%02.0f',rand*31+1),'',prm.dirs.nav);
for t=0:900:86400, tr=[tr;eceftogeod(navtostate(td,t,nav))]; end

% another matlab process -------------------------------------------------------
function execmd(cmd)
dos(['matlab -nojvm -automation -r ',cmd,' &']);
