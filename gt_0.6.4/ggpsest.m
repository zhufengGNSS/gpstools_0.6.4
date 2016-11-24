function h=ggpsest(varargin)
%-------------------------------------------------------------------------------
% [system] : GpsTools
% [module] : gpsestd gui interface
% [func]   : gpsestd gui interface
% [argin]  : none
% [argout] : h = figure handle
% [note]   :
% [version]: $Revision: 16 $ $Date: 2008-12-12 15:49:30 +0900 (金, 12 12 2008) $
%            Copyright(c) 2004-2008 by T.Takasu, all rights reserved
% [history]: 04/11/15  0.1  new
%            05/06/28  0.2  support ver.0.5.5
%            06/06/29  0.3  support ver.0.6.3
%            07/01/06  0.4  fix bug on replacing keyword in directories (gt_0.6.3p4)
%            08/11/21  0.5  support igs30s,cod5s,igscod5s for sat clock (gt_0.6.4)
%                           support itrf2005,igs05 for rcv position
%                           support vmf1 for tropspheric mapping function
%                           support shapiro-delay for relativistic correction
%                           support keyword in file path for output directory
%                           support phase-adjust to repair clock-jump
%                           support uigetfile bug of matlab 7.4
%-------------------------------------------------------------------------------
if nargin<1, h=GpsEstGui; else feval(varargin{:}); end

% gpsest gui interface ---------------------------------------------------------
function h=GpsEstGui
sel1={{'OFF','ON'},0:1};
sel2={{'RINEX','RINEX(3H)','RINEX(1H)','RINEX(15min)'},...
      {'rinex','rinex3','rinex1','rinexh'}};
sel3={{'ZD/PPP','DD'},{'zerodiff','dbldiff'}};
sel4={{'Combined','RINEX','RINEX(3H)','RINEX(1H)','RINEX(15min)'},...
      {'brdc','rinex','rinex3','rinex1','rinexh'}};
sel5={{'OFF','1pass(F)','2pass(FB)','3pass(FBF)'},0:3};
sel6={{'White Noise','Gauss-Marcov'},0:1};
sel7={{'Static','Kinematic'},0:1};
sel8={{'LC','L1','L2','L1/L2'},{'LC','L1','L2','L1L2'}};
prm1={
'','Estimation Start/End Time (GPST)', ''
't','',                                 [2000,1,1,0,0,0]
't','',                                 [2000,1,1,0,0,0]
'n','Estimation Interval (sec)',        300
's','Est Unit/Overlap Time (hr)',       [24,0]
'p','Exec Obs Data Editor',             sel1
'p','Exec Parameter Estimator',         sel5
'p','Estimation Strategy',              sel3
'b','Est./Fixed Parameters',            [mfilename,' cbEst1']
'b','Baselines',                        [mfilename,' cbEst2']
'p','Satellite Clock Model',            sel6
'p','Receiver Clock Model',             sel6
'p','Receiver Position Model',          sel7
'p','Clock Reference',                  {{''},{}}
};
prm2={
'p','Raw Observation Data',             sel2
'p','Navigation Messages',              sel4
'b','Obs Data Editor Settings',         [mfilename,' cbObsSetting']
'b','Estimation/Meas. Model',           [mfilename,' cbEst5']
'b','Measurement Noises',               [mfilename,' cbEstB']
'b','A Priori Std. Deviations',         [mfilename,' cbEst6']
'b','Process Noises',                   [mfilename,' cbEst7']
'b','Satellite Orbit Model',            [mfilename,' cbEst8']
'b','Quality Control Settings',         [mfilename,' cbEst9']
'b','Data Directories/Files',           [mfilename,' cbEst4']
'','',''
};
prm3={
' ','Clean Observation Data Directory','';
'd','',''
' ','Output Directory','';
'd','',''
};
h=gut('newfig','gtest','Obs Data Editor/Parameter Estimator',[600,400,0,-50],'');
if isempty(h), return, end
set(gcbo,'userdata',h);
h1=gut('newbtn','',[ 3,378,73,20],'Satellites...',[mfilename,' cbSat']);
h2=gut('newbtn','',[80,378,73,20],'Receivers...',  [mfilename,' cbRcv']);
set(h1,'fontsize',8);
set(h2,'fontsize',8);
gut('newlist','sats',[3,3,75,373],{},{},2,'');
gut('newlist','rcvs',[80,3,75,373],{},{},2,'');
gut('newprms','prm1',[165,80,210,23],prm1);
gut('newprms','prm2',[384,126,210,23],prm2);
gut('newprms','prm3',[384,80,210,17],prm3);
set(gut('newtext','pmsg',[164,50,428,20],'',2),'foregroundcolor',[0 0 0.8]);
gut('newprog','prog',[164,34,428,15],[0,0.8,0.8]);
gut('newbtnh', 'btns',[162,4,433,22],...
    {'File','Plot','Log','Execute','Close'},...
    {[mfilename,' cbFile'],[mfilename,' cbPlot'],[mfilename,' cbLog'],...
     [mfilename,' cbExe'],[mfilename,' cbExit']});
gut('newcmenu','file',[162,136],...
    {'&Files...','&Download...','-','&Execute Batch...','-','&Load Setting...',...
     '&Save Setting...'},{},...
     {[mfilename,' cbFiles'],[mfilename,' cbDown'],'',[mfilename,' cbExecBatch'],...
      '',[mfilename,' cbLoad'],[mfilename,' cbSave']});
gut('newcmenu','plot',[250,186],...
    {'&Observation Data','-','&Satellite Orbit...','Satellite &Clock...',...
     'Receiver &Position...','Receiver C&lock...','&Tropospheric Parameters...',...
     '&Residuals/Statistics...','-','Satellite Tra&ck/Station Pos...'},{},...
    {[mfilename,' cbPlotObs' ],'',[mfilename,' cbPlotEph' ],[mfilename,' cbPlotClks'],...
     [mfilename,' cbPlotPos'],[mfilename,' cbPlotClkr'],[mfilename,' cbPlotTrop'],...
     [mfilename,' cbPlotStats'],'',[mfilename,' cbPlotTrack']});
set(gut('geth','prm1_8'),'Callback',[mfilename,' cbUpdateEna']);
LoadSetting;

% callback on menus ------------------------------------------------------------
function cbFile, gut('showmenu','file');
function cbDown, execmd('gpsdown_');
function cbPlot, gut('showmenu','plot');

function cbFiles
p=GetSetting; gpsfiles; gpsfiles('cbUpdate',p.dirs.est);

function cbPlotObs
p=GetSetting; prm=loadprm('prm_plotobs','prm_plotobs_def');
prm.tstart=p.tstart; prm.tend=p.tend; prm.tint=p.tint; prm.tunit=p.tunit;
prm.sats=p.sats; prm.rcvs=p.rcvs; prm.src=p.obs.src;
if isempty(p.dirs.obs), prm.dirs.obs=p.dirs.est; else prm.dirs.obs=p.dirs.obs; end
if isempty(p.dirs.obc), prm.dirs.obc=p.dirs.est; else prm.dirs.obc=p.dirs.obc; end
if isempty(p.dirs.nav), prm.dirs.nav=p.dirs.est; else prm.dirs.nav=p.dirs.nav; end
prm.dirs.obs=strrep(prm.dirs.obs,'%O',p.dirs.est);
prm.dirs.obc=strrep(prm.dirs.obc,'%O',p.dirs.est);
prm.dirs.nav=strrep(prm.dirs.nav,'%O',p.dirs.est);
prm.navsrc=p.src.nav;
prm.type={'Data Availability',prm.type{:}}; prm.plotf=3;
if ~isempty(prm.rcvs), prm.rcv=prm.rcvs{1}; else prm.rcv={}; end
if ~isempty(prm.sats), prm.sat=prm.sats{1}; else prm.sat={}; end
plotobs('prm',prm);

function cbPlotEph
p=GetSetting; prm=loadprm('prm_ploteph','prm_ploteph_def');
prm.tstart=p.tstart; prm.tend=p.tend; prm.tint=p.tint; prm.tunit=p.tunit;
prm.sats=p.sats; prm.dirs.est=p.dirs.est;
if p.backward, prm.fb='ephfb'; else prm.fb='ephf'; end
if ~isempty(prm.sats), prm.sat=prm.sats{1}; else prm.sat={}; end
ploteph('prm',prm);

function cbPlotPos
p=GetSetting; prm=loadprm('prm_plotpos','prm_plotpos_def');
prm.tstart=p.tstart; prm.tend=p.tend; prm.tint=p.tint; prm.tunit=p.tunit;
prm.rcvs=p.rcvs; prm.dirs.est=p.dirs.est;
if p.backward, prm.fb='posfb'; else prm.fb='posf'; end
if ~isempty(prm.rcvs), prm.rcv=prm.rcvs{1}; else prm.rcv={}; end
plotpos('prm',prm);

function cbPlotClks
p=GetSetting; prm=loadprm('prm_plotclk','prm_plotclk_def');
prm.tstart=p.tstart; prm.tend=p.tend; prm.tint=p.tint; prm.tunit=p.tunit;
prm.sr='sat'; prm.sats=p.sats; prm.dirs.est=p.dirs.est;
if p.backward, prm.fb='clkfb'; else prm.fb='clkf'; end
if ~isempty(prm.sats), prm.sat=prm.sats{1}; else prm.sat={}; end
plotclk('prm',prm);

function cbPlotClkr
p=GetSetting; prm=loadprm('prm_plotclk','prm_plotclk_def');
prm.tstart=p.tstart; prm.tend=p.tend; prm.tint=p.tint; prm.tunit=p.tunit;
prm.sr='rcv'; prm.rcvs=p.rcvs; prm.dirs.est=p.dirs.est;
if p.backward, prm.fb='clkfb'; else prm.fb='clkf'; end
if ~isempty(prm.rcvs), prm.rcv=prm.rcvs{1}; else prm.rcv={}; end
plotclk('prm',prm);

function cbPlotTrop
p=GetSetting; prm=loadprm('prm_plottrop','prm_plottrop_def');
prm.tstart=p.tstart; prm.tend=p.tend; prm.tint=p.tint; prm.tunit=p.tunit;
prm.rcvs=p.rcvs; prm.dirs.est=p.dirs.est;
if p.backward, prm.fb='zpdfb'; else prm.fb='zpdf'; end
if ~isempty(prm.rcvs), prm.rcv=prm.rcvs{1}; else prm.rcv={}; end
plottrop('prm',prm);

function cbPlotStats
p=GetSetting; prm=loadprm('prm_plotstats','prm_plotstats_def');
prm.tstart=p.tstart; prm.tend=p.tend; prm.tint=p.tint; prm.tunit=p.tunit;
prm.sats=p.sats; prm.rcvs=p.rcvs(1); prm.sat='ALL'; prm.rcv='ALL'; prm.type='res';
prm.dirs.est=p.dirs.est;
if isempty(p.dirs.nav), prm.dirs.nav=p.dirs.est; else prm.dirs.nav=p.dirs.nav; end
prm.srcnav=p.src.nav;
if p.backward, prm.fb='resb'; else prm.fb='resf'; end
plotstats('prm',prm);

function cbPlotTrack
p=GetSetting; prm=loadprm('prm_gpstrack','prm_gpstrack_def');
prm.tstart=p.tstart; prm.tend=p.tend; prm.tint=max(p.tint,300); prm.elmin=p.elmin;
prm.sats=p.sats; prm.rcvs=p.rcvs; prm.sat='ALL'; prm.rcv='ALL';
if isempty(p.dirs.nav), prm.dirs.nav=p.dirs.est; else prm.dirs.nav=p.dirs.nav; end
prm.dirs.nav=strrep(prm.dirs.nav,'%O',p.dirs.est);
prm.src.nav=p.src.nav;
gpstrack('prm',prm);

% set setting ------------------------------------------------------------------
function SetSetting(prm)
gut('setstring','sats',prm.sats);
gut('setstring','rcvs',prm.rcvs);
gut('setlist','prm1_14',{'',prm.rcvs{:}},{});
pass=prm.gpsest+prm.backward+prm.iteration;
gut('setprms','prm1','',prm.tstart,prm.tend,prm.tint,[prm.tunit,prm.tover],...
    prm.obsedit,pass,prm.omodel,'','',prm.sclkmodel,prm.rclkmodel,...
    prm.rposmodel,prm.clkref);
gut('setprms','prm2',prm.obs.src,prm.src.nav);
gut('setprms','prm3','',prm.dirs.obc,'',prm.dirs.est);
gut('setudata','gtest',prm);
cbUpdateEna;

% get setting ------------------------------------------------------------------
function prm=GetSetting
prm=gut('getudata','gtest'); if isempty(prm), prm=prm_gpsest_def; end
prm.sats=gut('getstring','sats');
prm.rcvs=gut('getstring','rcvs');
[q,prm.tstart,prm.tend,prm.tint,tunit,prm.obsedit,pass,prm.omodel,q,q,...
 prm.sclkmodel,prm.rclkmodel,prm.rposmodel,prm.clkref]=gut('getprms','prm1');
prm.gpsest=pass>0; prm.backward=pass>1; prm.iteration=pass>2;
prm.tunit=tunit(1); prm.tover=tunit(2);
[prm.obs.src,prm.src.nav]=gut('getprms','prm2');
[q,prm.dirs.obc,q,prm.dirs.est]=gut('getprms','prm3');
gut('setudata','gtest',prm);

% update enable/disable ---------------------------------------------------------
function cbUpdateEna
prm=GetSetting;
e1=any(strcmp(prm.omodel,{'dbldiff','rcvdiff'}));
if e1, gut('setena','prm1_10'); else gut('setdis','prm1_10'); end
if ~isempty(prm.satest)
    e2=any([prm.satest{:,2}]==1); e3=any([prm.satest{:,4}]==1);
    if e2, gut('setena','prm2_8'); else gut('setdis','prm2_8'); end
    if e3, gut('setena','prm1_14'); else gut('setdis','prm1_14'); end
    if e3, gut('setena','prm1_11'); else gut('setdis','prm1_11'); end
end
if ~isempty(prm.rcvest)
    e4=any([prm.rcvest{:,2}]==1); e5=any([prm.rcvest{:,4}]==1);
    if e4, gut('setena','prm1_12'); else gut('setdis','prm1_12'); end
    if e5, gut('setena','prm1_13'); else gut('setdis','prm1_13'); end
end

% callback on satellite list ----------------------------------------------------
function cbSat
[sats,ok]=editlist('Satellites',gut('getstring','sats'));
if ok, gut('setlist','sats',sats,{}); end

% callback on receiver list -----------------------------------------------------
function cbRcv
prm=gut('getudata','gtest');
[rcvs,ok]=editlist('Receivers',gut('getstring','rcvs'),'rcv');
if ok
    gut('setlist','rcvs',rcvs,{});
    gut('setlist','prm1_14',{'',rcvs{:}},{});
    gut('setsel','prm1_14',prm.clkref);
end

% callback on load.../save... --------------------------------------------------
function cbLoad, settings(0);
function cbSave, settings(1);

function settings(save)
persistent path, prm=gut('getudata','gtest');
if isempty(path), path=fullfile(prm.dirs.est,'prm_*.mat'); end
pats={'*.mat','settings (prm_*.mat)';'*.m','settings (prm_*.m)'};
if save
    [file,p]=uiputfile(pats(1,:),'Save Settings',path);
    if file~=0, path=[p,file]; SaveSettingToFile(path); end
else
    [file,p]=uigetfile(pats,'Load Settings',path);
    if file~=0, path=[p,file]; LoadSettingFromFile(path); end
end

% callback on log... -----------------------------------------------------------
function cbLog
ShowLog(GetSetting);

% callback on execute... -------------------------------------------------------
function cbExe
if ~strcmp(gut('getstring','btns_4'),'Execute') % abort
    set(gut('geth','btns_4'),'enable','off','userdata',1);
    return
end
prm=SaveSetting;
hs=[gut('geth','pmsg'),gut('geth','prog'),gut('geth','btns_4')];
set(hs(3),'string','Abort','userdata',0);
f=gcf; gut('setdisall',f);
gut('setenah',hs(1)); gut('setenah',hs(3));
gmsg('seth',hs);
SaveHist(prm);
[stat,msg]=gpsestd(prm);
ms={'error : ','completed','aborted'};
gmsg('seth',[]);
set(hs(1),'string',[ms{stat+2},msg]);
set(hs(3),'string','Execute');
gut('setenaall',f);
figure(f); cbUpdateEna;

% callback on execute batch ----------------------------------------------------
function cbExecBatch
hs=[gut('geth','pmsg'),gut('geth','prog')];
prm=loadprm('prm_estbatch','prm_estbatch_def');
prm1={
't','Estimation Start Time (GPST)', prm.tstart
't','Estimation End Time (GPST)', prm.tend
'b','Receivers', [mfilename,' cbBatchRcv']
};
prm2={'d','', prm.dirs};
gut('newdlg','gtbatch','Execute Batch',[395,238]);
gut('newchk','sw1',[10,195,30,22],'',[mfilename,' cbBatchSw']);
gut('newchk','sw2',[10,160,30,22],'',[mfilename,' cbBatchSw']);
gut('newchk','sw3',[10,136,100,22],'Output Directory',[mfilename,' cbBatchSw']);
gut('newprms','prm1',[28,162,360,24,80],prm1);
gut('newprms','prm2',[119,139,270,24],prm2);
gut('newtext','',[10,113,300,20],'Processing Parameters');
gut('newlist','prms',[8,30,383,89],prm.prms,{},2,'');
h=gut('newbtnh', 'btns',[10,4,380,22],{'Add...','Delete','Log','Execute','Close'},...
      {[mfilename,' cbBatchAdd'],[mfilename,' cbBatchDel'],[mfilename,' cbBatchLog'],...
       [mfilename,' cbBatchExe'],[mfilename,' cbBatchClose']});
gut('setval','sw1',prm.sws(1));
gut('setval','sw2',prm.sws(2));
gut('setval','sw3',prm.sws(3));
gut('setudata','prm1_3',prm.rcvs);
gmsg('seth',[hs,h(4)]);
cbBatchSw;

function cbBatchSw
s1={'prm1_1_1','prm1_1_2','prm1_1_3','prm1_1_4','prm1_1_5','prm1_1_6',...
    'prm1_2_1','prm1_2_2','prm1_2_3','prm1_2_4','prm1_2_5','prm1_2_6'};
s2={'prm1_3'}; s3={'prm2_1','prm2_1b'};
if gut('getval','sw1'), gut('setena',s1); else gut('setdis',s1); end
if gut('getval','sw2'), gut('setena',s2); else gut('setdis',s2); end
if gut('getval','sw3'), gut('setena',s3); else gut('setdis',s3); end

function cbBatchRcv
rcvs=gut('getudata','prm1_3');
[rcvs,ok]=editlist('Receivers',rcvs,'rcv');
if ok, gut('setudata','prm1_3',rcvs); end

function cbBatchAdd
persistent path, if isempty(path), path=fullfile(pwd,'prm_*.mat'); end
[file,p]=uigetfile({'*.mat','settings (prm_*.mat)'},'Load Settings',path);
if file==0, return, end
path=[p,file]; prms=gut('getlist','prms'); gut('setlist','prms',{prms{:},path});

function cbBatchDel
prms=gut('getlist','prms');
for sel=gut('getsel','prms');, prms(strcmp(prms,sel{:}))=[]; end
gut('setlist','prms',prms);

function cbBatchLog
p=GetBatchPrms; if isempty(p.prms), return; end
prm=loadprm(p.prms{end},'prm_gpsest_def');
if p.sws(1), prm.tstart=p.tstart; prm.tend=p.tend; end
if p.sws(3), prm.dirs.est=p.dirs; end
ShowLog(prm);

function cbBatchExe
if ~strcmp(gut('getstring','btns_4'),'Execute') % abort
    set(gut('geth','btns_4'),'enable','off','userdata',1); return
end
f=gcf;
gut('setdisall',f);
set(gut('geth','btns_4'),'string','Abort','userdata',0);
gut('setena','btns_4');
p=GetBatchPrms; msg='';
for n=1:length(p.prms)
    if exist(p.prms{n})
        gut('setsel','prms',p.prms{n});
        prm=loadprm(p.prms{n},'prm_gpsest_def');
        if p.sws(1), prm.tstart=p.tstart; prm.tend=p.tend; end
        if p.sws(2), prm.rcvs=p.rcvs; end
        if p.sws(3), prm.dirs.est=p.dirs; end
        SaveHist(prm);
        [stat,msg]=gpsestd(prm); if stat~=0, break; end
    end
end
gut('setenaall',f); figure(f);
gut('setstring','btns_4','Execute');
ms={'error : ','completed','aborted'};
gmsg([ms{stat+2},msg]);

function cbBatchClose
prm=GetBatchPrms;
[path,file]=fileparts(which(mfilename));
try, save(fullfile(path,'settings','prm_estbatch.mat'),'prm'); catch ; end
closereq

function prm=GetBatchPrms
[prm.tstart,prm.tend]=gut('getprms','prm1');
prm.rcvs=gut('getudata','prm1_3');
prm.dirs=gut('getprms','prm2');
prm.sws=[gut('getval','sw1'),gut('getval','sw2'),gut('getval','sw3')];
prm.prms=gut('getlist','prms');

% callback on exit -------------------------------------------------------------
function cbExit
SaveSetting;
closereq

% obsedit setting dialog -------------------------------------------------------
function cbObsSetting
sel1={{'OFF','ON'},[0,1]};
sel2={{'OFF','Adjust Phase','Adjust Code'},0:2};
sel4={{'C1','P1','P1/C1'},{'C1','P1','P1C1'}};
prm1={
'n','L1 Frequency (GHz)',               1.57542
'n','L2 Frequency (GHz)',               1.22760
'n','Obs. Time Tolerance (sec)',        0.01
'p','Exclude Unhealthy Satellites',     sel1
'n','Min Elevation Angle (deg)',        0
'n','Max Arc Gap (sec)',                300
'n','Min Arc Length (sec)',             900
'n','Min Arc Points (points)',          10
'n','Noise of MW Code (m)',             0.5
'n','Noise of Ion-Free Code (m)',       1.5
'n','Slip Threshold f LG (m)',          0.03
'n','Slip Threshold of WL (cycle)',     0.5
'n','Slip Threshold of MP1 (m)',        0.5
'n','Slip Threshold of MP2 (m)',        0.5
};
prm2={
'n','Outlier Threshold of WL (sigma)',  4
'n','Ionos Smoothing Window (sec)',     0
'n','Moving Aver. Window 1 (pnts)',     60
'n','Moving Aver. Window 2 (pnts)',     60
'p','Slip Thres Elevation Weighting',   sel1
'p','Repair Steered Clock Jump',        sel2
'p','Repair Cycle Slip',                sel1
's','# of Pnts/Degrees for Poly. Fitt', [10,2]
'p','Separate Arc at Eclipse In/Out',   sel1
'p','Phase Multipath Correction',       sel1
's','MP Start/End/Offset (days/s)',     [0,0,0]
};
prm3={
' ','Code Multipath Profile',''
'f','',''
};
diss={'prm2_10','prm2_11_1_1','prm2_11_1_2','prm2_11_1_3','prm2_12',...
      '_prm2_10','_prm2_11','prm3_1','_prm3_1'}; 

h=gut('newdlg','gtedit','Observation Data Editor Settings',[520,335]);
gut('newprms','prm1',[10,8,245,23],prm1);
gut('newprms','prm2',[269,77,245,23],prm2);
gut('newprms','prm3',[269,31,245,17],prm3);
gut('setdis',diss);
prm=gut('getudata','gtest');
gut('setprms','prm1',prm.f1,prm.f2,prm.ttol,prm.exuhsat,prm.obs.elmin,...
    prm.obs.gapmax,prm.obs.arcmin,prm.obs.pntmin,prm.obs.sigmw,prm.obs.sigif,...
    prm.obs.dgmax,prm.obs.dwl,prm.obs.dmp1,prm.obs.dmp2)
gut('setprms','prm2',prm.obs.outl,prm.ionwind,prm.obs.wind(1),prm.obs.wind(2),...
    prm.obs.elwe,prm.obs.clkrep,prm.obs.reps,[prm.obs.npnt,prm.obs.nmax],...
    prm.obs.separc,prm.obs.srfilt,[prm.obs.srdays,prm.obs.sradj])
gut('setprms','prm3','',prm.mpc,'')
gut('newokcancelbtn','',[334,4,180,22]);
if ~gut('waitok'), return, end
prm.rcvs=gut('getstring','rcvs');
[prm.f1,prm.f2,prm.ttol,prm.exuhsat,prm.obs.elmin,prm.obs.gapmax,...
 prm.obs.arcmin,prm.obs.pntmin,prm.obs.sigmw,prm.obs.sigif,prm.obs.dgmax,...
 prm.obs.dwl,prm.obs.dmp1,prm.obs.dmp2]=gut('getprms','prm1');
[prm.obs.outl,prm.ionwind,prm.obs.wind(1),prm.obs.wind(2),prm.obs.elwe,...
 prm.obs.clkrep,prm.obs.reps,npnt,prm.obs.separc,prm.obs.srfilt,srdays...
]=gut('getprms','prm2');
prm.obs.npnt=npnt(1); prm.obs.nmax=npnt(2);
prm.obs.srdays=srdays(1:2); prm.obs.sradj=srdays(3);
[q,prm.mpc]=gut('getprms','prm3');
gut('setudata','gtest',prm);
close

% est/fixed parameters dialogs -------------------------------------------------
function cbEst1
sats=gut('getstring','sats'); sats={'','ALL',sats{:}};
rcvs=gut('getstring','rcvs'); rcvs={'','ALL',rcvs{:}};
sel1={'','<estimated>'};
sel2={'','<estimated>','IGS Final','IGS Rapid','IGS URapid','Est(forward)','Est(backward)','Est(smoothed)','ephs','IGS URapid(Pred)','CODE','EMR','ESA','GFZ','JPL','MIT','CODE Rapid'};
sel3={'','<estimated>','IGS Final','IGS Rapid','IGS URapid','Est(forward)','Est(backward)','Est(smoothed)','IGS URapid(Pred)','CODE','EMR','ESA','GFZ','JPL','MIT','IGS/CODE','CODE Rapid','IGR/CODER','IGS 30s','CODE 5s','IGS/CODE 5s'};
sel4={'','<estimated>','IGS Final','ITRF2000','ITRF97','GSIPOS','Est(forward)','Est(backward)','Est(smoothed)','IGS00','IGb00','IGS05','ITRF2005','Parameter File'};
sel5={'','<estimated>','IGS Final','IGS Rapid','IERS BulB','IERS BulA','IERS C04','Est(forward)','Est(backward)','Est(smoothed)','CODE','EMR','ESA','GFZ','JPL','MIT'};
sel6={'','IGS Final','IGS Rapid','IERS BulB','IERS C04','Est(forward)','Est(backward)','Est(smoothed)','CODE','EMR','ESA','GFZ','JPL','MIT'};
sel7={'','<estimated>','IGS ZPD','IGS ZPD Monthly','Est(forward)','Est(backward)','Est(smoothed)','JMA MSM Online'};
sel8={'','Est(forward)','Est(backward)','Est(smoothed)'};
val6={'','igs','igr','bulb','c04','erpf','erpb','erpfb','cod','emr','esa','gfz','jpl','mit'};
h=gut('newdlg','','Estimated/Fixed Parameters',[493,400]);
gut('newtexth','',[8,376,385,20],{'Satellite','Orbit','Clock','SRP Param.'},2);
gut('newtexth','',[8,214,480,20],{'Receiver','Position','Clock','Tropos. ZTD/Grad','Phase Bias'},2);
for n=1:6
    y=383-n*24;
    gut('newpopm',['satn_',num2str(n)],[  8,y,95,23],sats,{},'');
    gut('newpopm',['sato_',num2str(n)],[108,y,95,23],sel2,0:16,'');
    gut('newpopm',['satc_',num2str(n)],[203,y,95,23],sel3,0:20,'');
    gut('newpopm',['satr_',num2str(n)],[298,y,95,23],sel1,0:1,'');
end
for n=1:8
    y=221-n*24;
    gut('newpopm',['rcvn_',num2str(n)],[  8,y,95,23],rcvs,{},'');
    gut('newpopm',['rcvp_',num2str(n)],[108,y,95,23],sel4,0:13,'');
    gut('newpopm',['rcvc_',num2str(n)],[203,y,95,23],sel3(1:15),0:14,'');
    gut('newpopm',['rcvz_',num2str(n)],[298,y,95,23],sel7,0:7,'');
    gut('newpopm',['rcvb_',num2str(n)],[393,y,95,23],sel8,0:3,'');
end
gut('newtext','',[393,376,95,20],'ERP/Geocent',2);
gut('newtext','',[393,304,95,20],'Initial ERP ',2);
gut('newtext','',[393,256,95,20],'Input Unit Time(hr)',2);
gut('newpopm','erp',[393,359,95,23],sel5,0:15,'');
gut('newpopm','eco',[393,335,95,23],sel1,0:1,'');
gut('newpopm','erp0',[393,287,95,23],sel6,val6,'');
gut('newedit','tuinp',[393,239,95,23],'');
gut('newtext','',[8,2,100,20],'Day of Pos',1);
gut('newtext','',[190,2,100,20],'Interp Clock',1);
gut('newymd','tpos',[65,5,120,23],[],'');
gut('newpopm','clki',[248,3,50,23],{'OFF','ON'},0:1,'');
gut('newokcancelbtn','',[308,4,180,22]);
prm=gut('getudata','gtest');
for n=1:min(size(prm.satest,1),6)
    gut('setsel',['satn_',num2str(n)],prm.satest{n,1});
    gut('setsel',['sato_',num2str(n)],prm.satest{n,2});
    gut('setsel',['satc_',num2str(n)],prm.satest{n,4});
    gut('setsel',['satr_',num2str(n)],prm.satest{n,3});
end
for n=1:min(size(prm.rcvest,1),8)
    gut('setsel',['rcvn_',num2str(n)],prm.rcvest{n,1});
    gut('setsel',['rcvp_',num2str(n)],prm.rcvest{n,4});
    gut('setsel',['rcvc_',num2str(n)],prm.rcvest{n,2});
    gut('setsel',['rcvz_',num2str(n)],prm.rcvest{n,3});
    gut('setsel',['rcvb_',num2str(n)],prm.rcvest{n,5});
end
gut('setsel','erp',prm.est.erp);
gut('setsel','eco',prm.est.eco);
gut('setsel','erp0',prm.initerp);
gut('setymd','tpos',prm.tpos);
gut('setnum','tuinp',prm.tuinp);
gut('setsel','clki',prm.clkintp);
if ~gut('waitok'), return, end
prm.satest={}; m=0;
for n=1:6
    sat=gut('getsel',['satn_',num2str(n)]);
    if ~isempty(sat)
        m=m+1;
        prm.satest{m,1}=sat;
        prm.satest{m,2}=gut('getsel',['sato_',num2str(n)]);
        prm.satest{m,4}=gut('getsel',['satc_',num2str(n)]);
        prm.satest{m,3}=gut('getsel',['satr_',num2str(n)]);
    end
end
prm.rcvest={}; m=0;
for n=1:8
    rcv=gut('getsel',['rcvn_',num2str(n)]);
    if ~isempty(rcv)
        m=m+1;
        prm.rcvest{m,1}=rcv;
        prm.rcvest{m,4}=gut('getsel',['rcvp_',num2str(n)]);
        prm.rcvest{m,2}=gut('getsel',['rcvc_',num2str(n)]);
        prm.rcvest{m,3}=gut('getsel',['rcvz_',num2str(n)]);
        prm.rcvest{m,5}=gut('getsel',['rcvb_',num2str(n)]);
    end
end
prm.est.erp=gut('getsel','erp');
prm.est.eco=gut('getsel','eco');
prm.initerp=gut('getsel','erp0');
prm.tpos=gut('getymd','tpos');
prm.tuinp=gut('getnum','tuinp');
prm.clkintp=gut('getsel','clki');
gut('setudata','gtest',prm);
close;
cbUpdateEna;

% baseline for DD/TT dialog ---------------------------------------------------
function cbEst2
rcvs=gut('getstring','rcvs');
gut('newdlg','','Baselines for DD/TT',[380,290]);
gut('newokcancelbtn','',[195,4,180,22]);
prm=gut('getudata','gtest');
gut('newtext','',[10,265,165,20],'Receiver1    -    Receiver1',2);
gut('newtext','',[210,265,165,20],'Receiver1    -    Receiver2',2);
for n=1:10
    y=270-n*24;
    gut('newpopm', ['srcv_', num2str(n)],[10,y,80,23],{'',rcvs{:}},{},'');
    gut('newpopm', ['ercv_', num2str(n)],[95,y,80,23],{'','ALL',rcvs{:}},{},'');
    gut('newpopm', ['srcv_', num2str(n+10)],[210,y,80,23],{'',rcvs{:}},{},'');
    gut('newpopm', ['ercv_', num2str(n+10)],[295,y,80,23],{'','ALL',rcvs{:}},{},'');
end
for n=1:min(size(prm.baseline,1),20)
    gut('setsel', ['srcv_', num2str(n)],prm.baseline{n,1});
    gut('setsel', ['ercv_', num2str(n)],prm.baseline{n,2});
end
if ~gut('waitok'), return, end
prm.baseline={}; m=0;
for n=1:20
    srcv=gut('getsel',['srcv_', num2str(n)]);
    ercv=gut('getsel',['ercv_', num2str(n)]);
    if ~isempty(srcv)&~isempty(ercv)
        m=m+1; prm.baseline{m,1}=srcv; prm.baseline{m,2}=ercv;
    end
end
gut('setudata','gtest',prm);
close

% excluded sats/stas dialog ---------------------------------------------------
function cbEst3
sats=gut('getstring','sats');
rcvs=gut('getstring','rcvs');
[q,tstart,q,tspan]=gut('getprms','prm1');
list={'',sats{:},rcvs{:}};
gut('newdlg','','Excluded Satellites/Receivers',[492,338]);
gut('newtext','',[8,310,70,20],'Sats./Stas.',2);
gut('newtexth','',[82,310,405,20],{'Start Date/Time(GPST)','End Date/Time(GPST)'},2);
for n=1:12
    y=317-n*24;
    gut('newpopm', ['name_', num2str(n)],[10,y,70,23],{'',sats{:},rcvs{:}},{},'');
    gut('newedit', ['syear_',num2str(n)],[82,y+1,45,22],'');
    gut('newtable',['stime_',num2str(n)],[127,y+1,155,22],1,5);
    gut('newedit', ['eyear_',num2str(n)],[287,y+1,45,22],'');
    gut('newtable',['etime_',num2str(n)],[332,y+1,155,22],1,5);
end
h=gut('newbtn','',[10,3,80,22],'NANU...',[mfilename,' cbNanu']);
set(h,'userdata',{tstart,tspan});
gut('newokcancelbtn','',[307,4,180,22]);
prm=gut('getudata','gtest');
for n=1:min(size(prm.exclude),12)
    gut('setsel',  ['name_', num2str(n)],prm.exclude{n,1});
    gut('setnum',  ['syear_',num2str(n)],prm.exclude{n,2}(1));
    gut('settable',['stime_',num2str(n)],prm.exclude{n,2}(2:6));
    gut('setnum',  ['eyear_',num2str(n)],prm.exclude{n,3}(1));
    gut('settable',['etime_',num2str(n)],prm.exclude{n,3}(2:6));
end
if ~gut('waitok'), return, end
prm.exclude={}; m=0;
for n=1:12
    name=gut('getsel',['name_',num2str(n)]);
    if ~isempty(name)
        m=m+1;
        prm.exclude{m,1}=name;
        prm.exclude{m,2}=[gut('getnum',['syear_',num2str(n)]),...
                          gut('gettable',['stime_',num2str(n)])];
        prm.exclude{m,3}=[gut('getnum',['eyear_',num2str(n)]),...
                          gut('gettable',['etime_',num2str(n)])];
    end
end
gut('setudata','gtest',prm);
close

% apply scheduled mentenance form NANU ----------------------------------------
function cbNanu
persistent pp, if isempty(pp), pp='*.*'; end
t=get(gcbo,'userdata'); ts=caltomjd(t{1}); te=ts+t{2}/24; exc=[];
[file,path]=uigetfile(pp,'Load NANU'); if file==0, return, end
f=fopen([path,file],'rt'); if f==-1, return, end
while 1
    str=fgets(f); if ~isstr(str), break, end
    p=findstr(str,'UNUSABLE DAY');
    if ~isempty(p)
        s=sscanf(str(34:38),'PRN%d');
        t=sscanf(str(p+12:end),'%d/%d-%d/%d');
        if ~isempty(s)&length(t)>=3
            y=sscanf(str(1:4),'%d');
            if length(t)==3, t(3:4)=[t(1),t(3)]; end
            t1=daytime(y,t(1:2));
            t2=daytime(y,t(3:4));
            if caltomjd(t1)<te&ts<caltomjd(t2), exc=[exc;s,t1,t2]; end
        end
    end
end
fclose(f); pp=[path,'*.*'];
for n=1:size(exc,1)-1
    for m=n+1:size(exc,1)
        if all(exc(n,1:4)==exc(m,1:4)), exc(n,1)=0; break, end
    end
end
exc=exc(find(exc(:,1)),:);

for n=1:min(size(exc,1),12)
    gut('setsel',  ['name_', num2str(n)],sprintf('GPS%02d',exc(n,1)));
    gut('setnum',  ['syear_',num2str(n)],exc(n,2));
    gut('settable',['stime_',num2str(n)],exc(n,3:7));
    gut('setnum',  ['eyear_',num2str(n)],exc(n,8));
    gut('settable',['etime_',num2str(n)],exc(n,9:13));
end

function td=daytime(y,t)
td=mjdtocal(caltomjd([y,1,1])+t(1)-1,floor(t(2)/100)*3600+mod(t(2),100)*60);
td(6)=0;

% estimation/measurment model diag ---------------------------------------------
function cbEst5
rcvs=gut('getstring','rcvs');
sel1={{'No Tropos','Saastamoinen'},{'','trop_saast'}};
sel2={{'','COSZ','NMF','GMF','VMF1'},{'','mapf_cosz','mapf_nmf','mapf_gmf','mapf_vmf1'}};
sel3={{'OFF','ON'},[0,1]};
sel4={{'Standard Kalman'},{'filterekf'}};
sel6={{'IAU1976/1980','IERS1996'},{'','iers1996'}};
sel7={{'Standard Atmos.','GPT','JMA GSM Online','JMA MSM Online'},{'','gpt','gso','mso'}};
sel8={{'','Linear','Quadratic'},0:2};
sel9={{'Fixed Station','LEO Sattellite'},0:1};
sel10={{'Random-Walk','Gauss-Marcov'},0:1};
sel11={{'OFF','IGS Final','IGS Rapid','IGS-CODE'},{'','igs','igr','igscod'}};
sel12={{'OFF','Sat Clock','SatClk+Shapiro-Delay'},0:2};
prm1={
's','Min/Max Elevation Angle (deg)',    [10,90]
'p','Tropospheric A Priori Model',      sel1
'p','Tropospheric Mapping Function',    sel2
'p','Tropospheric Stochastic Model',    sel10
'p','Tropospheric Gradient Model',      sel8
'n','Meas. Noise Elev. Angle Factor',   1
's','Meas. Noise Factor in Eclipse',    [5,5]
's','Process Noise Factor in Eclipse',  [5,5]
'n','Max Post-Eclipse Maneuver (sec)',  3600
'n','Std. Dev. Factor of Initial States', 2
'n','Std. Dev. Factor of PreDet States', 0
'n','Max Count of Filter Iteration',     2
'n','Min Count of Obs for Estimation',  1
'p','Nutation Model',                   sel6
'p','Filter',                           sel4
'p','Clock Correction by Phase-Bias',   sel3
'p','QC/Align Clock to Reference',      sel11
};
prm2={
'p','Receiver Attitude Model',          sel9
' ','Measurement Corrections',          ''
'p','  Satellite Antenna Offset',       sel3
'p','  Receiver Antenna Offset',        sel3
'p','  Relativistic Effects',           sel12
'p','  Phase Windup',                   sel3
'p','  Phase Multipath',                sel3
' ','Station Displacement Corrections',  ''
'p','  Solid Earth Tides',              sel3
'p','  Ocean Loading',                  sel3
'p','  Pole Tides',                     sel3
'p','  Eliminate Permanent Deform.',    sel3
'p','ERP Variation Correction',         sel3
'p','Meteorological Parameters',        sel7
'n','Default MSL Pressure (hPa)',       0
's','Default Temp/Rel Humidity (C,%)',  [0,0]
};
gut('newdlg','gtest5','Estimation/Measurement Model',[550,398]);
gut('newprms','prm1',[15,5,260,23,95],prm1);
gut('newprms','prm2',[285,28,260,23,95],prm2);
gut('newokcancelbtn','',[366,4,180,22]);
prm=gut('getudata','gtest'); q='';
gut('setprms','prm1',[prm.elmin,prm.elmax],prm.trop,prm.mapf,prm.zpdmodel,...
    prm.trgmodel,prm.sigweight,prm.eclnsig,prm.eclnprn,prm.ecltime,prm.sdfact,...
    prm.sdfact2,prm.maxiter,prm.minobs,prm.nutmodel,prm.filter,prm.clkcorr,...
    prm.clkqc);
gut('setprms','prm2',prm.rattmodel,q,prm.corrf(1),prm.corrf(2),prm.corrf(3),...
    prm.corrf(4),prm.corrf(5),q,prm.sitedisp(1),prm.sitedisp(2),prm.sitedisp(3),...
    prm.sitedisp(4),prm.erpvar,prm.metsrc,prm.metprm(1),prm.metprm(2:3));
e1=~isempty(prm.rcvest)&any([prm.rcvest{:,3}]==1);
e3=~isempty(prm.rcvest)&any([prm.rcvest{:,2}]==1);
e2=~isempty(prm.satest)&any([prm.satest{:,4}]==1);
if e1, gut('setena',{'prm1_4','prm1_5'}); else gut('setdis',{'prm1_4','prm1_5'}); end
if e2, gut('setena',{'prm1_17'}); else gut('setdis',{'prm1_17'}); end
if e2|e3, gut('setena',{'prm1_16'}); else gut('setdis',{'prm1_16'}); end
if ~gut('waitok'), return, end
[el,prm.trop,prm.mapf,prm.zpdmodel,prm.trgmodel,prm.sigweight,prm.eclnsig,...
 prm.eclnprn,prm.ecltime,prm.sdfact,prm.sdfact2,prm.maxiter,prm.minobs,...
 prm.nutmodel,prm.filter,prm.clkcorr,prm.clkqc]=gut('getprms','prm1');
prm.elmin=el(1); prm.elmax=el(2); prm.clkalign=prm.clkqc;
[prm.rattmodel,q,prm.corrf(1),prm.corrf(2),prm.corrf(3),prm.corrf(4),...
 prm.corrf(5),q,prm.sitedisp(1),prm.sitedisp(2),prm.sitedisp(3),prm.sitedisp(4),...
 prm.erpvar,prm.metsrc,prm.metprm(1),prm.metprm(2:3)]=gut('getprms','prm2');
gut('setudata','gtest',prm);
close

% measurment noise dialog ------------------------------------------------------
function cbEstB
rcvs=gut('getstring','rcvs');
gut('newdlg','gtestb','Measurement Noises',[300,165]);
gut('newtext','',[15,138,270,23],'Measurement Noise (m)',2);
for n=1:5
    y=145-23*n;
    gut('newpopm',['rcvn_',num2str(n)],[12,y,70,23],{'','ALL',rcvs{:}},{},'');
    gut('newedit',['osig_',num2str(n)],[82,y+1,65,22],'');
    gut('newpopm',['rcvn_',num2str(n+5)],[158,y,70,23],{'','ALL',rcvs{:}},{},'');
    gut('newedit',['osig_',num2str(n+5)],[228,y+1,65,22],'');
end
gut('newokcancelbtn','',[113,4,180,22]);
prm=gut('getudata','gtest');
for n=1:min(size(prm.obssig),10)
    gut('setsel',['rcvn_',num2str(n)],prm.obssig{n,1});
    gut('setnum',['osig_',num2str(n)],prm.obssig{n,2});
end
if ~gut('waitok'), return, end
prm.obssig={};
for n=1:10
    rcv=gut('getsel',['rcvn_',num2str(n)]);
    if ~isempty(rcv)
        prm.obssig{n,1}=rcv;
        prm.obssig{n,2}=gut('getnum',['osig_',num2str(n)]);
    end
end
gut('setudata','gtest',prm);
close

% apriori std. deviations dialog -----------------------------------------------
function cbEst6
title='A Priori Standard Deviations';
sats=gut('getstring','sats'); sats={'','ALL',sats{:}};
rcvs=gut('getstring','rcvs'); rcvs={'','ALL',rcvs{:}};
prm=gut('getudata','gtest');
if length(prm.sig.eco)<2, prm.sig.eco(2)=5; end
prmopt=[prm.sig.erp,prm.sig.eco];
[prm.satsig,prm.rcvsig,prmopt]=...
    satrcvprmdlg(title,[550,400],sats,rcvs,prm.satsig,prm.rcvsig,prmopt,'',prm);
prm.sig.erp=prmopt(1:2); prm.sig.eco=prmopt(3:4);
gut('setudata','gtest',prm);

% process noise dialog ---------------------------------------------------------
function cbEst7
title='Process Noises';
sats=gut('getstring','sats'); sats={'','ALL',sats{:}};
rcvs=gut('getstring','rcvs'); rcvs={'','ALL',rcvs{:}};
prm=gut('getudata','gtest');
if length(prm.prn.eco)<2, prm.prn.eco(2)=1E-3; end
prmopt=[prm.prn.erp,prm.prn.eco];
[prm.satprn,prm.rcvprn,prmopt]=...
    satrcvprmdlg(title,[550,400],sats,rcvs,prm.satprn,prm.rcvprn,prmopt,...
                 'Units : 1/sqrt(sec)',prm);
prm.prn.erp=prmopt(1:2); prm.prn.eco=prmopt(3:4);
gut('setudata','gtest',prm);

% sat/sta setting dialog -------------------------------------------------------
function [prmsat,prmrcv,prmopt]=satrcvprmdlg(title,siz,sats,rcvs,prmsat,prmrcv,...
                                             prmopt,text,prm)
label1={'Pos(m)','Vel(m/s)','Srp1','Srp2','Srp3','Srp4-6','Clk(m)','dClk(m/h)',...
        'XpYp(rad)'};
label2={'Clk(m)','dClk(m/h)','ZTD(m)','dZTD(m/h)','Grad(m)','PosU(m)',...
        'PosE(m)','PosN(m)','PBias(m)'};
gut('newdlg','satrcvprm',title,siz);
gut('newtext','',[10,372,70,23],'Satellite',2);
gut('newtext','',[10,211,70,23],'Receiver',2);
gut('newtexth','',[80,372,465,23],label1,2);
gut('newtexth','',[80,211,465,23],label2,2);
gut('newtext','',[495,327,50,20],'Lod(sec)',2);
gut('newtext','',[495,292,50,20],'Gco',2);
gut('newtext','',[495,279,50,20],'(m/ppb)',2);
gut('newtext','',[88,2,100,20],text,2);
for n=1:6
    y=382-n*24;
    h1=gut('newpopm', ['sat_',num2str(n)],[10,y,70,23],sats,{},'');
    h2=gut('newtable',['san_',num2str(n)],[80,y+1,413.3,22],1,8);
    if any([prm.satest{:,:}]==1), e='on'; else e='off'; end
    set(h1,'enable',e); set(h2,'enable',e);
end
for n=1:8
    y=221-n*24;
    h1=gut('newpopm', ['rcv_',num2str(n)],[10,y,70,23],rcvs,{},'');
    h2=gut('newtable',['rcn_',num2str(n)],[80,y+1,465,22],1,9);
    if any([prm.rcvest{:,:}]==1), e='on'; else e='off'; end
    set(h1,'enable',e); set(h2,'enable',e);
end
h1=gut('newedit','erp1',[495,359,50,22],'');
h2=gut('newedit','erp2',[495,311,50,22],'');
if prm.est.erp==1, e='on'; else e='off'; end
set(h1,'enable',e); set(h2,'enable',e);
h1=gut('newedit','eco1',[495,263,50,22],'');
h2=gut('newedit','eco2',[495,239,50,22],'');
if prm.est.eco==1, e='on'; else e='off'; end
set(h1,'enable',e); set(h2,'enable',e);
gut('newokcancelbtn','',[367,4,180,22]);
for n=1:min(size(prmsat,1),6)
    gut('setsel',  ['sat_',num2str(n)],prmsat{n,1});
    gut('settable',['san_',num2str(n)],[prmsat{n,2:end}]);
end
for n=1:min(size(prmrcv,1),8)
    gut('setsel',  ['rcv_',num2str(n)],prmrcv{n,1});
    gut('settable',['rcn_',num2str(n)],[prmrcv{n,2:end}]);
end
gut('setnum','erp1',prmopt(1));
gut('setnum','erp2',prmopt(2));
gut('setnum','eco1',prmopt(3));
gut('setnum','eco2',prmopt(4));

if ~gut('waitok'), return, end
prmsat={}; prmrcv={};
for n=1:6
    sat=gut('getsel',['sat_',num2str(n)]);
    if ~isempty(sat)
        prmsat{n,1}=sat;
        sv=gut('gettable',['san_',num2str(n)]);
        for m=1:length(sv), prmsat{n,m+1}=sv(m); end
    end
end
for n=1:8
    rcv=gut('getsel',['rcv_',num2str(n)]);
    if ~isempty(rcv)
        prmrcv{n,1}=rcv;
        sv=gut('gettable',['rcn_',num2str(n)]);
        for m=1:length(sv), prmrcv{n,m+1}=sv(m); end
    end
end
prmopt(1)=gut('getnum','erp1');
prmopt(2)=gut('getnum','erp2');
prmopt(3)=gut('getnum','eco1');
prmopt(4)=gut('getnum','eco2');
close

% satellite orbit model dialog -------------------------------------------------
function cbEst8
sel1={{'OFF','Sun/Moon'},{'','grv_sunmoon'}};
sel2={{'OFF','Simple','ROCK4','GSPM.04','GSPM.04+A/C','CODE(D/Y/B)','CODE(D/Y/B/Z)','CODE(D/Y/B/Z/X)'},...
      {'','srp_simple','srp_rock4','srp_gspm','srp_gspmm','srp_code','srp_code2','srp_code3'}};
sel3={{'Cylindric','Penumbra/Umbra'},{'','penumbra'}};
sel4={{'OFF'},{''}};
sel5={{'OFF','Solid Earth Tide','Solid/Ocean Tide'},{'','solid','otide'}};
sel6={{'OFF','IERS'},{'','iers'}};
sel7={{'OFF','ON'},{'','dv_on'}};
prm1={
'n','Max Degrees of Geopotential',    8
'p','Solar/Planetary Potentials',     sel1
'p','Solar Radiation Pressure Model', sel2
'p','Eclipse Model',                  sel3
'p','Atmospheric Drag Model',         sel4
'p','Earth Tides',                    sel5
'p','Relativistic Effects',           sel6
'p','Thrust Force',                   sel7
'n','Integration Step (sec)',         30
};
gut('newdlg','gtest6','Satellite Orbit Model',[295,250]);
gut('newprms','prm1',[15,35,275,23,105],prm1);
gut('newokcancelbtn','',[112,4,180,22]);
prm=gut('getudata','gtest'); q='';
gut('setprms','prm1',prm.obt.g_nmax,prm.obt.p_plgrv,prm.obt.p_solarpr,...
    prm.obt.p_eclipse,prm.obt.p_atmos,prm.obt.p_tidal,prm.obt.p_relativ,...
    prm.obt.p_deltav,prm.obt.tstep);
if ~gut('waitok'), return, end
[prm.obt.g_nmax,prm.obt.p_plgrv,prm.obt.p_solarpr,prm.obt.p_eclipse,...
 prm.obt.p_atmos,prm.obt.p_tidal,prm.obt.p_relativ,prm.obt.p_deltav,...
 prm.obt.tstep]=gut('getprms','prm1');
gut('setudata','gtest',prm);
close

% quality control dialog -------------------------------------------------------
function cbEst9
sel1={{'OFF','ON'},0:1};
prm1={
'p','Re-Estimate without Excluded',        sel1
'n','Outlier Thres. of Prefit Res.(sigma)', 4
'n','Outlier Thres. of Postfit Res.(sigma)',4
'n','Arc Max Outlier Rate',                0
'n','Arc Max Postfit Residual (m)',        0
'n','Satellite Max Outlier Rate',          0
'n','Satellite Max Postfit Residual (m)',  0
'n','Receiver Max Outlier Rate',           0
'n','Receiver Max Postfit Residual (m)',   0
'n','Sat. Max Outage of Valid Obs (sec)',  3600
'n','Rcv. Max Outage of Valid Obs (sec)',  10800
'n','Satellite Clock Max Std Dev. (ns)',   0
'n','Receiver Clock Max Std Dev. (ns)',    0
'n','Receiver Position Max Std Dev. (m)',  0
'p','Output Statistics/Residuals', sel1
'p','Debug Error Stop/Output Debug Info', sel1
};
gut('newdlg','gtest9','Quality Control Settings',[297,390]);
gut('newprms','prm1',[15,33,275,22,90],prm1);
gut('newokcancelbtn','',[112,4,180,22]);
prm=gut('getudata','gtest'); q='';
gut('setprms','prm1',prm.reests,prm.outlp,prm.outlf,prm.varc.pout,prm.varc.rmsf,...
    prm.vsat.pout,prm.vsat.rmsf,prm.vrcv.pout,prm.vrcv.rmsf,prm.maxsatout,...
    prm.maxrcvout,prm.maxsatcsig,prm.maxrcvcsig,prm.maxrcvpsig,prm.statout,prm.dbout);
if ~gut('waitok'), return, end
[prm.reests,prm.outlp,prm.outlf,prm.varc.pout,prm.varc.rmsf,prm.vsat.pout,...
 prm.vsat.rmsf,prm.vrcv.pout,prm.vrcv.rmsf,prm.maxsatout,prm.maxrcvout,...
 prm.maxsatcsig,prm.maxrcvcsig,prm.maxrcvpsig,prm.statout,prm.dbout]=gut('getprms','prm1');
gut('setudata','gtest',prm);
close

% data directories/files dialog ------------------------------------------------
function cbEst4
label1={'Raw Observation Data','Navigation Messages','Satellite Ephemerides',...
        'Satellite/Receiver Clocks','Receiver Positions','Earth Rotation Parameters',...
        'Tropos/Meteo Parameters','Input Estimated Sat Params.',...
        'Input Estimated Rcv Params.'};
label2={'Satellite Antenna PCV','Receiver Antenna PCV','Receiver Position',...
        'Ocean Loading Parameter','Excluded Sats/Rcvs','Phase Multipath Profile'};
prm1={
'd','','';'d','','';'d','','';'d','','';'d','','';'d','','';'d','','';'d','','';...
'd','',''};
prm2={'h','','';'h','','';'h','','';'h','','';'h','','';'f','',''};
gut('newdlg','gtest4','Data Directories/Files',[480,360]);
gut('newtextv','',[15,166,160,21],label1);
gut('newtextv','',[15,32,160,21],label2);
gut('newprms','prm1',[155,164,320,21],prm1);
gut('newprms','prm2',[155,30,320,21],prm2);
info=sprintf('Matlab %s, %s, %s',version,lower(computer),version('-java'));
set(gut('newtext','',[15,10,240,10],info),'foregroundcolor',[.5 .5 .5],'fontsize',7);
set(gut('newbtn','',[260,6,16,18],'?',[mfilename,' cbPathHelp']),'fontsize',8);
gut('newokcancelbtn','',[295,4,180,22]);
prm=gut('getudata','gtest'); q='';
gut('setprms','prm1',prm.dirs.obs,prm.dirs.nav,prm.dirs.eph,prm.dirs.clk,...
    prm.dirs.pos,prm.dirs.erp,prm.dirs.trop,prm.dirs.inp,prm.dirs.inpr);
gut('setprms','prm2',prm.satpcv,prm.pcv,prm.rcvposf,prm.oload,prm.exsatrcv,...
    prm.mpp);
if ~gut('waitok'), return, end
[prm.dirs.obs,prm.dirs.nav,prm.dirs.eph,prm.dirs.clk,prm.dirs.pos,...
 prm.dirs.erp,prm.dirs.trop,prm.dirs.inp,prm.dirs.inpr]=gut('getprms','prm1');
[prm.satpcv,prm.pcv,prm.rcvposf,prm.oload,prm.exsatrcv,prm.mpp]=gut('getprms','prm2');
gut('setudata','gtest',prm);
close

function cbPathHelp
str1={'%S','%s','%G','%Y','%y','%m','%d','%n','%D','%W','%P','%O'};
str2={': Sat/Rcv Name (UPPER Case)',': Sat/Rcv Name (lower Case)',...
      ': Sat/Rcv Name (Tail 4-Chars)',': Year (yyyy)',': Year (yy)',...
      ': Month (mm)',': Day of Month (dd)',': Day of Year (ddd)',...
      ': Day of Week (d)',': GPS Week No (wwww)',': Install Directory',...
      ': Output Directory'};
gut('newdlg','','Keyword Replacements',[220,218]);
gut('newtextv','',[20,20,30,16],str1);
gut('newtextv','',[50,20,200,16],str2);
gut('newokbtn','',[145,4,70,23]);
gut('waitok');
close;

% load setting -----------------------------------------------------------------
function LoadSetting
SetSetting(loadprm('prm_gpsest','prm_gpsest_def'));

% save setting -----------------------------------------------------------------
function prm=SaveSetting
[dirs,f]=fileparts(which(mfilename));
prm=SaveSettingToFile(fullfile(dirs,'settings','prm_gpsest.mat'));

% load setting from file -------------------------------------------------------
function LoadSettingFromFile(file)
if exist(file), SetSetting(loadprm(file,'prm_gpsest_def')); end

% save setting to file ---------------------------------------------------------
function prm=SaveSettingToFile(file)
prm=GetSetting;
try, save(file,'prm'); catch gut('errdlg',lasterr); end

% save processing history ------------------------------------------------------
function SaveHist(prm)
[dirs,f]=fileparts(which(mfilename));
epoch=datevec(now);
file=sprintf('prm_%04d%02d%02d%02d%02d%02.0f.mat',epoch);
try, save(fullfile(dirs,'history',file),'prm'); catch ; end

% show log ---------------------------------------------------------------------
function ShowLog(prm)
if ~isempty(prm.rcvs), name=prm.rcvs{1}; else name=''; end
log={}; dirs=gfilepath(prm.dirs.est,'',prm.tstart,name);
for f=dir(fullfile(dirs,'gpsestd*.log'))'
    log={log{:},fullfile(dirs,f.name)};
end
gut('newviewer','','Parameter Estimator Log ',[600,450],'file',unique(log));

% another matlab process -------------------------------------------------------
function execmd(cmd)
v=sscanf(version,'%f');
if v(1)<=7.3
    dos(['matlab -nojvm -automation -r ',cmd,' &']);
else
    dos(['matlab -automation -r ',cmd,' &']); % support v.7.4 uicontrol bug
end