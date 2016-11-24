function h=ploteph(varargin)
%-------------------------------------------------------------------------------
% [system] : GpsTools
% [module] : plot satellite orbits
% [func]   : plot satellite orbits
% [argin]  : 'prm',prm   = parameters
% [argout] : h = figurea handle
% [note]   :
% [version]: $Revision: 20 $ $Date: 06/07/18 14:03 $
%            Copyright(c) 2004-2006 by T.Takasu, all rights reserved
% [history]: 04/12/22  0.1  new
%            05/01/19  0.2  add erp plot function
%            05/04/11  0.3  add eco plot function
%-------------------------------------------------------------------------------
if nargin>0&strncmp(varargin{1},'cb',2), feval(varargin{:});
else h=PlotMain(varargin{:}); end

% main -------------------------------------------------------------------------
function h=PlotMain(varargin)
h=gut('newfigg','','Satellite Orbit',[600,400,0,-68],[mfilename,' cbClose']);
if isempty(h), return, end
prm=loadprm('prm_ploteph','prm_ploteph_def');
if nargin>=2&strcmp(varargin{1},'prm'), prm=varargin{2}; end
cb1=[mfilename,' cbPlot1'];
cb2=[mfilename,' cbPlot2'];
cb3=[mfilename,' cbPlot3'];
cb4=[mfilename,' cbPlot4'];
gut('newmenu','data','&Data ',...
    {'&Read...','&New Window','-','Out&put...','&Export Plot...','-','&Options...'},{},...
    {[mfilename,' cbRead'],[mfilename,' cbNew'],'',[mfilename,' cbOut'],...
     [mfilename,' cbExport'],'',[mfilename,' cbOpt']});
gut('newmenu','plot','&Plot',...
    {'Position Error &3D','Position Error &R/A/C','Position Error &X/Y/Z','-',...
     'Position Error by &Satellites','Position Error by &Dates','-',...
     '&SRP Parameters','&Earth Rotation Parameters','&Geocenter Offset','-',...
     'Summary...','Output SRP &Parameters...'},...
    {'pos3d','posrac','posxyz','','possat','posdat','','srp','erp','eco','','',''},...
    {cb1,cb1,cb1,'',cb2,cb2,'',cb2,cb2,cb2,'',cb3,cb4});
gut('newmenu','sats','&Satellite',prm.sats,prm.sats,[mfilename,' cbSat']);
gut('chkmenu','plot',prm.type);
gut('chkmenu','sats',prm.sat);
set(gut('newbtnh','sbtn',[0,0,36,16],{'<','>'},[mfilename,' cbSwBtn']),...
    'fontsize',8);
set(gcf,'resizefcn',[mfilename,' cbResize']);
LocatePos;
data.eph={}; data.sig={}; data.err={};
data.srp={}; data.sigsrp={};
data.erp={}; data.sigerp={};
data.eco={}; data.sigeco={};
data.prm=prm; set(h,'userdata',data);
if length(varargin)>0, ReadData; UpdatePlot; end

% callback on new --------------------------------------------------------------
function cbNew
data=get(gcf,'userdata');
feval(mfilename,'prm',data.prm);

% callback on close ------------------------------------------------------------
function cbClose
SaveSetting;
closereq

% callback on resize -----------------------------------------------------------
function cbResize, LocatePos

% callback on menu read --------------------------------------------------------
function cbRead
e1={'Estimated(Forward)','Estimated(Backward)','Estimated(Smoothed)'};
e2={'ephf','ephb','ephfb'};
r1={'IGS Final','IGS Rapid','IGS URapid','IGS URapid(Pred)','EPHG','EPHS','COD','EMR','ESA','GFZ','JPL','MIT','NGS','SIO','CODE Rapid','Broardcast'};
r2={'igs','igr','igu','igp','ephg','ephs','cod','emr','esa','gfz','jpl','mit','ngs','sio','codr','brdc'};
data=get(gcf,'userdata'); prm=data.prm;
prm1={
't','Start Time (GPST)',[]
't','End Time (GPST)',  []
'n','Time Interval (sec)',0
'n','Processing Unit Time (hr)',24
'b','Satellites',[mfilename,' cbSelSat']
'p','Satellite Orbit',{{e1{:},r1{:}},{e2{:},r2{:}}}
'p','Reference Satellite Orbit',{r1,r2}
'p','Interpolate Reference Orbit',{{'OFF','ON'},0:1}
};
prm2={
' ','Satellite Orbit Directory',''
'd','',prm.dirs.est
' ','Reference Satellite Orbit Directory',''
'd','',prm.dirs.eph
' ','Earth Rotation Parameters Directory',''
'd','',prm.dirs.erp
};
gut('newdlg','','Read Data',[314,352]);
gut('newprms','prm1',[12,146,295,25,118],prm1);
gut('newprms','prm2',[12,33,295,18],prm2);
gut('newokcancelbtn','',[127,4,180,23]);
gut('setprms','prm1',prm.tstart,prm.tend,prm.tint,prm.tunit,'',prm.fb,prm.ref,prm.interp);
gut('setudata','prm1_5',prm.sats);
if ~gut('waitok'), return, end
[prm.tstart,prm.tend,prm.tint,prm.tunit,q,prm.fb,prm.ref,prm.interp]=gut('getprms','prm1');
[q,prm.dirs.est,q,prm.dirs.eph,q,prm.dirs.erp]=gut('getprms','prm2');
prm.sats=gut('getudata','prm1_5');
if isempty(prm.sats), prm.sat=''; else prm.sat=prm.sats{1}; end
close
data.prm=prm; set(gcf,'userdata',data);
SaveSetting; ReadData; UpdatePlot;

% callback on menu select satellite --------------------------------------------
function cbSelSat
[sats,ok]=editlist('Satellites',get(gcbo,'userdata'));
if ok, set(gcbo,'userdata',sats); end

% callback on menu output results ----------------------------------------------
function cbOut
data=get(gcf','userdata'); prm=data.prm;
i=find(strcmp(prm.sat,prm.sats));
if isempty(data.eph)|isempty(i), return, end
s1={};
s1=adds(s1,'%%      SATELLITE POSITION : %s : %s-%s GPST',prm.sats{i},tstr(prm.td,prm.time(1)),tstr(prm.td,prm.time(end)));
s1=adds(s1,'%%  date     time(sec)       x(m)           y(m)           z(m)       sdx(m)  sdy(m)  sdz(m)');
s1=adds(s1,'%%------------------------------------------------------------------------------------------');
s2=cell(1,length(prm.time));
for n=1:length(prm.time)
    d=mjdtocal(prm.td,prm.time(n)); t=mod(prm.time(n),86400);
    s2{n}=sprintf('%04d,%02d,%02d,%9.2f, %14.4f,%14.4f,%14.4f, %7.4f,%7.4f,%7.4f\n',...
                  d(1:3),t,data.eph{i}(n,1:3),data.sig{i}(n,1:3));
end
Viewer({s1{:},s2{:}},['Satellite Orbit : ',prm.dirs.est,' (',prm.fb,')']);

% callback on menu export plot -------------------------------------------------
function cbExport
persistent file dirs
if isempty(file), file='ploteph'; end, if isempty(dirs), dirs=''; end
f=gcf; data=get(f,'userdata'); prm=data.prm;
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
if isempty(strmatch(prm.type,'pos'))
    sat=prm.sat;
    for n=1:length(prm.rcvs)
        prm.sat=prm.sats{n}; data.prm=prm; set(f,'userdata',data)
        figure(f), UpdatePlot;
        print(format{:},opts{:},sprintf('%s%02d',fs,n));
    end
    prm.sat=sat; data.prm=prm; set(f,'userdata',data); UpdatePlot;
else
    print(format{:},opts{:},fs);
end

% callback on menu options -----------------------------------------------------
function cbOpt
sel1={{'OFF','ON'},0:1};
sel2={{'Dot','Line','Line and Dot'},1:3};
data=get(gcf,'userdata'); prm=data.prm;
prm1={
'p','Show Mean and RMS Error',        sel1
'n','Show Est. Std Dev (sigma,0:OFF)', 0
'p','Plot Style',                     sel2
's','Plot Line Width / Marker Size',  [0.1,10]
's','Plot Axis Range (m)',            [-0.2,0.2]
};
gut('newdlg','','Options',[314,165]);
gut('newprms','prm1',[12,34,294,25,100],prm1);
gut('newokcancelbtn','',[126,4,180,23]);
gut('setprms','prm1',prm.showf(1),prm.nsig,prm.ptype,prm.psize,prm.range);
if ~gut('waitok'), return, end
[prm.showf(1),prm.nsig,prm.ptype,prm.psize,range]=gut('getprms','prm1');
if range(1)<range(2), prm.range=range; else prm.range=[-0.2,0.2]; end
close
data.prm=prm; set(gcf,'userdata',data);
SaveSetting; UpdatePlot;

% callback on menu plot --------------------------------------------------------
function cbPlot1
data=get(gcf,'userdata'); prm=data.prm;
types={'possat','posdat','srp','erp','eco'};
gut('unchkmenu','plot',types);
gut('togglechk',gcbo);
prm.type=gut('getchk','plot');
data.prm=prm; set(gcf,'userdata',data); UpdatePlot;

function cbPlot2
data=get(gcf,'userdata'); prm=data.prm;
prm.type={get(gcbo,'userdata')};
gut('unchkmenu','plot');
gut('chkmenu','plot',prm.type);
data.prm=prm; set(gcf,'userdata',data); UpdatePlot;

function cbPlot3, ErrSummary
function cbPlot4, SrpSummary

% callback on menu satellite ---------------------------------------------------
function cbSat
data=get(gcf,'userdata'); prm=data.prm;
prm.sat=get(gcbo,'userdata');
gut('unchkmenu','sats');
gut('chkmenu','sats',prm.sat);
data.prm=prm; set(gcf,'userdata',data);
UpdatePlot

% callback on push swithch button ----------------------------------------------
function cbSwBtn
data=get(gcf,'userdata'); prm=data.prm;
j=find(strcmp(prm.sats,prm.sat));
switch get(gcbo,'tag')
case 'sbtn_1', prm.sat=prm.sats{max(j-1,1)};
case 'sbtn_2', prm.sat=prm.sats{min(j+1,length(prm.sats))};
end
gut('unchkmenu','sats'); gut('chkmenu','sats',prm.sat);
data.prm=prm; set(gcf,'userdata',data);
UpdatePlot

% locate position --------------------------------------------------------------
function LocatePos
p=get(gcf,'position'); m=15;
gut('setpos','sbtn_1',[p(3)-m-32,p(4)-16]);
gut('setpos','sbtn_2',[p(3)-m-16,p(4)-16]);

% save settingsrs --------------------------------------------------------------
function SaveSetting
data=get(gcf,'userdata'); prm=data.prm;
[path,file]=fileparts(which(mfilename));
save(fullfile(path,'settings','prm_ploteph.mat'),'prm');

% update plot ------------------------------------------------------------------
function UpdatePlot
data=get(gcf,'userdata'); prm=data.prm;
if isempty(data.eph), return, end
if     any(strcmp(prm.type,'erp')), PlotErp
elseif any(strcmp(prm.type,'srp')), PlotSrp
elseif any(strcmp(prm.type,'eco')), PlotEco
elseif any(strcmp(prm.type,'possat')), PlotEstErrSat
elseif any(strcmp(prm.type,'posdat')), PlotEstErrDate
else PlotEstErr, end

% plot estimation error --------------------------------------------------------
function PlotEstErr
data=get(gcf,'userdata'); prm=data.prm;
clf
i=find(strcmp(prm.sat,prm.sats));
if isempty(data.err)|isempty(data.sig)|isempty(i), return, end
err=data.err{i}; sig=data.sig{i};
if isempty(err), return, end
margin=[0.08,0.03,0.04,0.012];
for n=1:length(prm.type)
    if n<length(prm.type), topts='nolabel'; else topts=''; end
    ggt('subplotv','est',n,length(prm.type),'move',[1,1,1,1],'link',[1,0,1,0],...
        'taxis',prm.td,'topts',topts,'margin',margin,...
        'xlim',[prm.time(1),prm.time(end)]/3600);
    switch prm.type{n}
    case 'pos3d'
        i=find(~isnan(err(:,1)));
        ylim=yrange(prm.range,err(i,1));
        plotval(prm.time,err(:,1),prm.ptype,prm.psize);
        plotsig(prm.time,0,sig(:,1),prm.nsig,'r:');
        set(gca,'ylim',ylim-ylim(1));
        if prm.showf(1)
            s=sprintf('REF: %s  MEAN: %7.4fm RMS: %7.4fm',RefName(prm.ref),...
                      meann(err(:,1)),rmsn(err(:,1)));
            ggt('mtext',['est_',num2str(n)],s,1);
        end
        ylabel('Position Error 3D (m)')
    case 'posrac'
        i=find(~isnan(err(:,5)));
        set(gca,'ylim',yrange(prm.range,err(i,5:7)));
        h=plotval(prm.time,err(:,5:7),prm.ptype,prm.psize);
        plotsig(prm.time,0,sig(:,5:7),prm.nsig,':');
        if ~isempty(h)
            legend(h,{'Radial','Along-Track','Cross-Track'});
        end
        if prm.showf(1)
            s=sprintf(['REF: %s\nMEAN R:%7.4fm A:%7.4fm C:%7.4fm\n',...
                       'RMS R:%7.4fm A:%7.4fm C:%7.4fm'],...
                      RefName(prm.ref),meann(err(:,5:7)),rmsn(err(:,5:7)));
            ggt('mtext',['est_',num2str(n)],s,3);
        end
        ylabel('Position Error R/A/C (m)')
    case 'posxyz'
        i=find(~isnan(err(:,2)));
        set(gca,'ylim',yrange(prm.range,err(i,2:4)));
        h=plotval(prm.time,err(:,2:4),prm.ptype,prm.psize);
        plotsig(prm.time,0,sig(:,2:4),prm.nsig,':');
        if ~isempty(h)
            legend(h,{'X (m)','Y (m)','Z (m)'});
        end
        if prm.showf(1)
            s=sprintf(['REF: %s\nMEAN X:%7.4fm Y:%7.4fm Z:%7.4fm\n',...
                       'RMS X:%7.4fm Y:%7.4fm Z:%7.4fm'],...
                      RefName(prm.ref),meann(err(:,2:4)),rmsn(err(:,2:4)));
            ggt('mtext',['est_',num2str(n)],s,3);
        end
        ylabel('Position Error X/Y/Z (m)')
    end
    if n==1
        ti=[tstr(prm.td,prm.time(1)),'-',tstr(prm.td,prm.time(end))];
        title(['Satellite Orbit ',prm.sat,' : ',ti]);
    end
end

% plot estimation error by satellites ------------------------------------------
function PlotEstErrSat
data=get(gcf,'userdata'); prm=data.prm;
clf
pos=[0.08,0.12,0.89,0.81];
for n=1:length(prm.sats)
    err(n,:)=rmsn(data.err{n}(:,[1,5:7]));
end
err=[err;meann(err)];
h=ggt('barplot',err,{prm.sats{:},'Average'},'ylim',prm.range-prm.range(1),...
      'position',pos);
legend(h,{'3D','Radial','Along-Track','Cross-Track'})
ylabel('Position RMS Error (m)');
ti=['Satellite Orbit : ',tstr(prm.td,prm.time(1)),'-',tstr(prm.td,prm.time(end)),...
    ' (REF: ',RefName(prm.ref),')'];
title(ti);

% plot estimatin error by date -------------------------------------------------
function PlotEstErrDate
data=get(gcf,'userdata'); prm=data.prm;
clf
pos=[0.08,0.12,0.89,0.81];
tu=prm.tunit*3600; ts=(floor(prm.time(1)/tu):floor(prm.time(end)/tu))*tu;
err=[]; label={};
for n=1:length(ts)
    i=find(ts(n)<=prm.time&prm.time<=ts(n)+tu-prm.tint); e=[];
    for m=1:length(prm.sats)
        if ~isempty(data.err{m}), e=[e;rmsn(data.err{m}(i,[1,5:7]))]; end
    end
    err=[err;meann(e)];
    ep=mjdtocal(prm.td,ts(n)); label={label{:},sprintf('%d/%d',ep(2:3))};
end
err=[err;meann(err)]; label={label{:},'Average'};
h=ggt('barplot',err,label,'ylim',prm.range-prm.range(1),'position',pos);
legend(h,{'3D','Radial','Along-Track','Cross-Track'})
ylabel('Position RMS Error (m)');
ti=['Satellite Orbit : ',tstr(prm.td,prm.time(1)),'-',tstr(prm.td,prm.time(end)),...
    ' (REF: ',RefName(prm.ref),')'];
title(ti);

% plot erp --------------------------------------------------------------------
function PlotErp
label={'Xp (")','Yp (")','UT1-UTC (sec)'};
data=get(gcf,'userdata'); prm=data.prm;
clf
erp=data.erp; sig=data.sigerp;
tr=(floor((prm.time(1)-43200)/86400):floor((prm.time(end)+43200)/86400))*86400+43200;
for n=1:length(tr)
    erpr(n,:)=readerp(prm.td+tr(n)/86400,prm.dirs.erp,'igs');
    erpr(n,1:2)=erpr(n,1:2)*180/pi*3600;
end
ti=['Earth Rotation Parameters : ',tstr(prm.td,prm.time(1)),'-',tstr(prm.td,prm.time(end))];
if isempty(erp), return, end
margin=[0.08,0.03,0.04,0.012];
for n=1:3
    if n<3, topts='nolabel'; else topts=''; end
    ggt('subplotv','est',n,3,'move',[1,1,1,1],'link',[1,0,1,0],...
        'taxis',prm.td,'topts',topts,'margin',margin,...
        'xlim',[prm.time(1),prm.time(end)]/3600);
    i=find(~isnan(erp(:,n)));
    if ~isempty(i), set(gca,'ylim',yrange([-inf,inf],erp(i,n))); end
    h=plotval(prm.time,erp(:,n),prm.ptype,prm.psize);
    plotsig(prm.time,erp(:,n),sig(:,n),prm.nsig,'r:');
    ylabel(label{n})
    h=plot(tr/3600,erpr(:,n),'m.');
    if n==1, legend(h,'IGS Final'); title(ti); end
end

% plot geocenter offset --------------------------------------------------------------------
function PlotEco
label={'dx (m)','dy (m)','dz (m)','scale (ppb)'};
data=get(gcf,'userdata'); prm=data.prm;
clf
eco=data.eco; sig=data.sigeco;
ti=['Geocenter Offset : ',tstr(prm.td,prm.time(1)),'-',tstr(prm.td,prm.time(end))];
if isempty(eco), return, end
margin=[0.08,0.03,0.04,0.012];
for n=1:4
    if n<4, topts='nolabel'; else topts=''; end
    ggt('subplotv','est',n,4,'move',[1,1,1,1],'link',[1,0,1,0],...
        'taxis',prm.td,'topts',topts,'margin',margin,...
        'xlim',[prm.time(1),prm.time(end)]/3600);
    i=find(~isnan(eco(:,n)));
    if ~isempty(i)
        if n==4, set(gca,'ylim',yrange([-inf,inf],eco(i,n)));
        else set(gca,'ylim',[-0.01,0.01]); end
    end
    h=plotval(prm.time,eco(:,n),prm.ptype,prm.psize);
    plotsig(prm.time,eco(:,n),sig(:,n),prm.nsig,'r:');
    ylabel(label{n})
    if n==1, title(ti); end
end

% plot srp parameters ----------------------------------------------------------
function PlotSrp
data=get(gcf,'userdata'); prm=data.prm;
clf
i=find(strcmp(prm.sat,prm.sats));
if isempty(data.srp)|isempty(data.sig)|isempty(i), return, end
srp=data.srp{i}; sig=data.sigsrp{i};
ti=['SRP Parameters ',prm.sat,' : ',tstr(prm.td,prm.time(1)),'-',tstr(prm.td,prm.time(end))];
margin=[0.08,0.03,0.04,0.012];
ns=max(find(~isnan(srp(1,:))));
for n=1:ns
    if n<ns, topts='nolabel'; else topts=''; end
    ggt('subplotv','est',n,ns,'move',[1,1,1,1],'link',[1,0,1,0],...
        'taxis',prm.td,'topts',topts,'margin',margin,...
        'xlim',[prm.time(1),prm.time(end)]/3600);
    i=find(~isnan(srp(:,n)));
    if ~isempty(i), set(gca,'ylim',yrange([-inf,inf],srp(i,n))); end
    h=plotval(prm.time,srp(:,n),prm.ptype,prm.psize);
    plotsig(prm.time,srp(:,n),sig(:,n),prm.nsig,'r:');
    ylabel(sprintf('srp%d',n))
    if n==1, title(ti); end
end

% summary of satellite position error ------------------------------------------
function ErrSummary
data=get(gcf','userdata'); prm=data.prm;
if isempty(data.err), return, end
s={};
s=adds(s,'                SATELLITE ORBIT : POSITION ERROR (wrt %s)',RefName(prm.ref));
s=adds(s,'');
s=AddSummary(s,data,prm,prm.time(1),prm.time(end));
s=adds(s,'');
s=adds(s,' -----------------------------------------------------------------------------');
tu=prm.tunit*3600;
for ts=(floor(prm.time(1)/tu):floor(prm.time(end)/tu))*tu
    s=adds(s,'');
    s=AddSummary(s,data,prm,ts,ts+tu-prm.tint);
end
if isempty(prm.files), ti=[prm.dirs.est,' ',RefName(prm.fb)]; else ti=prm.files{1}; end
Viewer(s,ti);

function s=AddSummary(s,data,prm,ts,te)
s=adds(s,' TIME     : %s-%s GPST',tstr(prm.td,ts),tstr(prm.td,te));
s=adds(s,'                       rms error (m)                  mean error (m) (ecef)   ');
s=adds(s,'              3d     radial along-trk cross-trk    x        y        z        ');
s=adds(s,' -----------------------------------------------------------------------------');
i=find(ts<=prm.time&prm.time<=te);
for n=1:length(prm.sats)
    if ~isempty(data.err{n})
        e(n,:)=[rmsn(data.err{n}(i,[1,5:7])),meann(data.err{n}(i,2:4))];
    else
        e(n,:)=repmat(nan,1,7);
    end
    s=adds(s,' %-7s : %7.4f %8.4f %8.4f %8.4f  %8.4f %8.4f %8.4f',prm.sats{n},e(n,:));
end
s=adds(s,' -----------------------------------------------------------------------------');
s=adds(s,' %-7s : %7.4f %8.4f %8.4f %8.4f  %8.4f %8.4f %8.4f','average',meann(e));
s=adds(s,' %-7s : %7.4f %8.4f %8.4f %8.4f  %8.4f %8.4f %8.4f','median',mediann(e));

% summary srp parameters -------------------------------------------------------
function SrpSummary
data=get(gcf','userdata'); prm=data.prm;
if isempty(data.srp), return, end;
s={};
s=adds(s,'%%');
s=adds(s,'%%                     SRP PARAMETERS : ESTIMATED VALUES');
s=adds(s,'%%               (%s-%s GPST)',tstr(prm.td,prm.time(1)),tstr(prm.td,prm.time(end)));
s=adds(s,'%%');
s=adds(s,'function prms=prm_srpprms');
s=adds(s,'prms={');
form='''%s'',%9.4f'; for n=1:11, form=[form,',%7.4f']; end
for n=1:6, form=[form,',%7.1E']; end
tu=prm.tunit*3600;
for n=1:length(prm.sats)
    srp=[];
    for t=(floor(prm.time(1)/tu):floor(prm.time(end)/tu))*tu
        i=find(t<=prm.time&prm.time<t+tu);
        if isempty(data.srp{n}), srp=[srp;nan];
        elseif prm.fb=='ephf', srp=[srp;data.srp{n}(i(end),:)];
        else srp=[srp;data.srp{n}(i(1),:)]; end
    end
    sigsrp=stdn(srp);
    prnsrp=sqrt(rmsn(srp(1:end-1,:)-srp(2:end,:)).^2/tu); % 1/sqrt(sec)
    sigsrp(isnan(sigsrp))=0; prnsrp(isnan(prnsrp))=0;
    s=adds(s,form,prm.sats{n},meann(srp),sigsrp,prnsrp);
end
s=adds(s,'};');
if isempty(prm.files), ti=[prm.dirs.est,' ',RefName(prm.fb)]; else ti=prm.files{1}; end
Viewer(s,ti);

% read estimation data ---------------------------------------------------------
function ReadData
data=get(gcf,'userdata'); prm=data.prm;
[prm.td,ts]=caltomjd(prm.tstart); [tn,te]=caltomjd(prm.tend);
time=ts:prm.tint:te+(tn-prm.td)*86400; if isempty(time), return, end
prm.time=time;
data.eph={}; data.err={}; data.sig={}; data.srp={}; data.sigsrp={};
gut('newmsgbar','msgbar','',[480,50],[0,0.7,1]);
gut('setmsgbar','msgbar','reading satellite orbit',0);
if any(strcmp(prm.fb,{'ephf','ephb','ephfb'}))
    [data.eph,data.sig,data.srp,data.sigsrp]=...
        ReadEphD(prm.td,prm.time,prm.sats,prm.dirs.est,prm.fb,prm.tunit);
else
    [ephs,sigs]=readeph(prm.td,prm.time,prm.sats,prm.dirs.est,prm.fb);
    for n=1:length(prm.sats)
       data.eph{n}=ephs(:,:,n);
       data.sig{n}=sigs(:,:,n);
    end
end
ti=['Satellite Orbit : ',prm.dirs.est,' (',prm.fb,')'];

% read erp
if any(strcmp(prm.fb,{'ephf','ephb','ephfb'}))
    gut('setmsgbar','msgbar','reading erp',0.5);
    [t,xs,vs]=readest(prm.td,prm.time,'erp','',prm.dirs.est,prm.fb(4:end),prm.tunit);
    if ~isempty(xs)
        [tt,i,j]=intersect(prm.time,t);
        data.erp=repmat(nan,length(prm.time),3); data.sigerp=data.erp;
        data.erp(i,:)=xs(j,:); data.sigerp(i,:)=sqrt(vs(j,:));
    end
end
% read references and compute errors
gut('setmsgbar','msgbar','reading reference satellite orbit',0.75);
[data.err,data.sig]=ErrEst(prm.td,prm.time,prm.sats,prm.dirs,data.eph,data.sig,...
                           prm.ref,prm.tunit,prm.interp);

gut('closemsgbar','msgbar');
data.prm=prm; set(gcf,'name',ti,'userdata',data);

% update menu
delete(gut('geth','sats'));
gut('newmenu','sats','&Satellite',prm.sats,prm.sats,[mfilename,' cbSat']);
gut('chkmenu','sats',prm.sat);

% ephemeris estimation errors --------------------------------------------------
function [err,sig]=ErrEst(td,time,sats,dirs,eph,sigs,ref,tunit,interp)
err={}; sig={};
if interp
    ephr=readeph(td,time,sats,dirs.eph,ref,tunit,'interp');
else
    ephr=readeph(td,time,sats,dirs.eph,ref,tunit); % reference ephemeris (ecef)
end
for n=1:length(sats)
   err{n}=repmat(nan,length(time),7); sig{n}=err{n};
   if ~isempty(eph{n})&~isempty(ephr)
       err{n}(:,2:4)=eph{n}(:,1:3)-ephr(:,1:3,n);         % poserr(ecef)
       sig{n}(:,2:4)=sigs{n}(:,1:3);                      % possig(ecef)
       
       for m=1:length(time)
           err{n}(m,1)=norm(err{n}(m,2:4));               % poserr(3d)
           sig{n}(m,1)=norm(sig{n}(m,2:4));               % possig(3d)
           
           E=eceftosatf(eph{n}(m,:)');
           err{n}(m,5:7)=err{n}(m,2:4)*E';                % poserr(rac)
           sig{n}(m,5:7)=sqrt(diag(E*diag(sig{n}(m,2:4).^2)*E'))'; % possig(rac)
       end
   end
end

% read estimated ephemeris -----------------------------------------------------
function [eph,sig,srp,sigsrp]=ReadEphD(td,time,sats,dirs,ephsrc,tunit)
ephs=repmat(nan,[length(time),6,length(sats)]); sigs=ephs; prm=[];
eph={}; sig={}; srp={}; sigsrp={};
for n=1:length(sats), eph{n}=[]; sig{n}=[]; srp{n}=[]; sigsrp{n}=[]; end
tu=tunit*3600;
for ts=(floor(time(1)/tu):floor(time(end)/tu))*tu
    i=find(ts<=time&time<ts+tu);
    for n=1:length(sats)
        [t,xs,vs,p]=readest(td,time(i),'eph',sats{n},dirs,ephsrc(4:end),tunit);
        [t,j,k]=intersect(time(i),t);
        if ~isempty(j)
            if isempty(srp{n})
                srp{n}=repmat(nan,length(time),size(xs,2)-6); sigsrp{n}=srp{n};
            end
            ephs(i(j),:,n)=xs(k,1:6);
            sigs(i(j),:,n)=sqrt(vs(k,1:6));
            srp{n}(i(j),:)=xs(k,7:end);
            sigsrp{n}(i(j),:)=sqrt(vs(k,7:end));
        end
        if ~isempty(p), prm=p; end
    end
    if ~isempty(prm)
        [ephs(i,:,:),sigs(i,:,:)]=...
            ephtoecef(td,time(i),ephs(i,:,:),sigs(i,:,:),prm,sats,dirs,ephsrc(4:end),tunit);
    end
end
for n=1:length(sats)
    eph{n}=ephs(:,:,n);
    sig{n}=sigs(:,:,n);
end

% reference name ---------------------------------------------------------------
function name=RefName(ref)
s={'Estimated(Forward)','Estimated(Backward)','Estimated(Smoothed)',...
   'IGS Final','IGS Rapid','IGS URapid','IGS URapid(Pred)','EPHG','EPHS','COD','EMR','ESA','GFZ','JPL','MIT','NGS','SIO'};
c={'ephf','ephb','ephfb','igs','igr','igu','igp','ephg','ephs','cod','emr','esa','gfz','jpl','mit','ngs','sio'};
i=find(strcmp(ref,c));
if isempty(i), name=''; else name=s{i}; end

% show viewer ------------------------------------------------------------------
function Viewer(str,ti)
gut('newviewer','',['Summary of Satellite Orbit : ',ti],[600,419,0,-20],'str',str);

% plot value -------------------------------------------------------------------
function h=plotval(t,x,ptype,psize)
if ptype==1|ptype==3, h=plot(t/3600,x,'.','markersize',psize(2)); end
if ptype==2|ptype==3, h=plot(t/3600,x,'-','linewidth',psize(1)); end

% plot sigma -------------------------------------------------------------------
function plotsig(t,x,sig,nsig,pat)
if nsig>0, plot(t/3600,x-nsig*sig,pat), plot(t/3600,x+nsig*sig,pat), end

% range ------------------------------------------------------------------------
function ylim=yrange(range,x)
ylim=range;
if isinf(ylim(1)), ylim(1)=min(min(x)); end
if isinf(ylim(2)), ylim(2)=max(max(x)); end
if any(isnan(ylim)), ylim=[0,1]; end

% mean/rms without nan ---------------------------------------------------------
function m=meann(x)
for n=1:size(x,2)
    i=find(~isnan(x(:,n)));
    if isempty(i), m(1,n)=nan; else m(1,n)=mean(x(i,n),1); end
end
function m=rmsn(x)
for n=1:size(x,2)
    i=find(~isnan(x(:,n)));
    if isempty(i), m(1,n)=nan; else m(1,n)=sqrt(mean(x(i,n).^2,1)); end
end
function m=stdn(x)
for n=1:size(x,2)
    i=find(~isnan(x(:,n)));
    if isempty(i), m(1,n)=nan; else m(1,n)=std(x(i,n),1); end
end
function m=mediann(x)
for n=1:size(x,2)
    i=find(~isnan(x(:,n)));
    if isempty(i), m(1,n)=nan; else m(1,n)=sqrt(median(x(i,n).^2,1)); end
end

% date/time string -------------------------------------------------------------
function str=tstr(td,ts)
str=sprintf('%04d/%02d/%02d %02d:%02d:%02.0f',mjdtocal(td,ts));

% add string -------------------------------------------------------------------
function s=adds(s,varargin), s={s{:},sprintf(varargin{:})};
