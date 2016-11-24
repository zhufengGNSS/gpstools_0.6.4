function h=plotbcp(varargin)
%-------------------------------------------------------------------------------
% [system] : GpsTools
% [module] : plot phase bias estimation
% [func]   : plot phase bias estimation
% [argin]  : none
% [argout] : h = figure handle
% [note]   :
% [version]: $Revision: 12 $ $Date: 2008-11-25 10:02:15 +0900 (火, 25 11 2008) $
%            Copyright(c) 2004-2006 by T.Takasu, all rights reserved
% [history]: 06/03/14  0.1  new
%-------------------------------------------------------------------------------
if nargin>0&strncmp(varargin{1},'cb',2), feval(varargin{:})
else h=PlotMain(varargin{:}); end

% main -------------------------------------------------------------------------
function h=PlotMain(varargin)
h=gut('newfigg','','Phase Bias Parameters',[600,400,0,-68],[mfilename,' cbClose']);
if isempty(h), return, end
prm=loadprm('prm_plotbcp','prm_plotbcp_def');
if nargin>=2&strcmp(varargin{1},'prm'), prm=varargin{2}; end
cb1=[mfilename,' cbPlot1'];
gut('newmenu','data','&Data',...
    {'&Read...','Read &File...','&New Window','-','&Options...'},{},...
    {[mfilename,' cbRead'],[mfilename,' cbFile'],[mfilename,' cbNew'],'',...
     [mfilename,' cbOpt']});
gut('newmenu','plot','&Plot',{'&Phase Bias'},{'bcp'},{cb1});
gut('newmenu','rcvs','&Receiver',prm.rcvs,prm.rcvs,[mfilename,' cbRcv']);
gut('chkmenu','plot',prm.type);
gut('chkmenu','rcvs',prm.rcv);
set(gcf,'resizefcn',[mfilename,' cbResize']);
data.prm=prm; data.bcp={}; data.sig={}; sata.sats={}; set(h,'userdata',data);
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
function cbResize

% callback on menu read --------------------------------------------------------
function cbRead
sel1={{'Estimated(Forward)','Estimated(Backward)','Estimated(Smoothed)'},...
      {'bcpf','bcpb','bcpfb'}};
data=get(gcf,'userdata'); prm=data.prm;
prm1={
't','Start Time (GPST)',[]
't','End Time (GPST)',  []
'n','Time Interval (sec)', 300
'n','Processing Unit Time (hr)',  24
'b','Receivers',[mfilename,' cbSelRcv']
'p','Phase Bias Parameters',sel1
};
prm2={
' ','Phase Bias Parameters Directory',''
'd','',prm.dirs.est
};
gut('newdlg','','Read Data',[314,248]);
gut('newprms','prm1',[12,82,295,27,118],prm1);
gut('newprms','prm2',[12,38,295,18],prm2);
gut('newokcancelbtn','',[127,4,180,23]);
gut('setprms','prm1',prm.tstart,prm.tend,prm.tint,prm.tunit,'',prm.fb);
gut('setudata','prm1_5',prm.rcvs);
if ~gut('waitok'), return, end
[prm.tstart,prm.tend,prm.tint,prm.tunit,q,prm.fb]=gut('getprms','prm1');
[q,prm.dirs.est]=gut('getprms','prm2');
prm.rcvs=gut('getudata','prm1_5');
if ~isempty(prm.rcvs), prm.rcv=prm.rcvs{1}; end
close
data.prm=prm; set(gcf,'userdata',data);
SaveSetting; ReadData; UpdatePlot;

% callback on menu select receiver --------------------------------------------
function cbSelRcv
[rcvs,ok]=editlist('Receivers',get(gcbo,'userdata'),'rcv');
if ok, set(gcbo,'userdata',rcvs); end

% callback on menu read file --------------------------------------------------
function cbFile
persistent path, if isempty(path), path=''; end
data=get(gcf,'userdata'); prm=data.prm;
[f,p]=uigetfile(fullfile(path,'bcp*.mat'),'Read File'); if f==0, return, end
path=p; prm.files={[p,f]};
data.prm=prm; set(gcf,'userdata',data)
ReadData, UpdatePlot

% callback on menu options -----------------------------------------------------
function cbOpt
sel1={{'OFF','ON'},0:1};
sel2={{'Dot','Line','Line and Dot'},1:3};
data=get(gcf,'userdata'); prm=data.prm;
prm1={
'n','Show Est. Std Dev. (sigma,0:OFF)', 0
'p','Plot Line Type',                  sel2
's','Plot Line Width / Marker Size',   [0.5,5]
's','Plot Axis Range (m)',             [-0.5,0.5]
};
q='';
gut('newdlg','','Options',[320,134]);
gut('newprms','prm1',[12,34,300,23,120],prm1);
gut('newokcancelbtn','',[133,4,180,23]);
gut('setprms','prm1',prm.nsig,prm.ptype,prm.psize,prm.range)
if ~gut('waitok'), return, end
[prm.nsig,prm.ptype,prm.psize,range]=gut('getprms','prm1');
if range(1)<range(2), prm.range=range; else prm.range=[-0.5,0.5]; end
close
data.prm=prm; set(gcf,'userdata',data)
SaveSetting; UpdatePlot;

% callback on menu receiver ----------------------------------------------------
function cbRcv
data=get(gcf,'userdata'); prm=data.prm;
prm.rcv=get(gcbo,'userdata');
gut('unchkmenu','rcvs');
gut('chkmenu','rcvs',prm.rcv);
data.prm=prm; set(gcf,'userdata',data);
UpdatePlot;

% save settingsrs --------------------------------------------------------------
function SaveSetting
data=get(gcf,'userdata'); prm=data.prm;
[path,file]=fileparts(which(mfilename));
save(fullfile(path,'settings','prm_plotbcp.mat'),'prm');

% update plot ------------------------------------------------------------------
function UpdatePlot
PlotBcp;

% plot phase bias parameters ---------------------------------------------------
function PlotBcp
data=get(gcf,'userdata'); prm=data.prm;
if isempty(data.bcp), return, end
i=find(strcmp(prm.rcv,prm.rcvs));
bcp=data.bcp{i}; sig=data.sig{i};
clf
margin=[0.08,0.03,0.04,0.012];
ggt('subplotv','est',1,1,'move',[1,1,1,1],'link',[1,0,1,0],'taxis',prm.td,...
    'margin',margin,'xlim',[prm.time(1),prm.time(end)]/3600);
if ~isempty(bcp)
    plotest(prm.time,bcp,sig,prm.nsig,prm.ptype,prm.psize);
end
ylabel('Phase Bias (m)');
ti=[tstr(prm.td,prm.time(1)),'-',tstr(prm.td,prm.time(end))];
title(['Phase Bias Parameters ',prm.rcv,' : ',ti]);

% read bcp data ----------------------------------------------------------------
function ReadData
data=get(gcf,'userdata'); prm=data.prm;
data.bcp={}; data.sig={}; data.sats={};
[prm.td,ts]=caltomjd(prm.tstart); [tn,te]=caltomjd(prm.tend);
prm.time=ts:prm.tint:te+(tn-prm.td)*86400; if isempty(prm.time), return, end
h=gcf; set(h,'pointer','watch');
for n=1:length(prm.rcvs)
    data.bcp{n}=[]; data.sig{n}=[];
    [t,xs,vs,p,arc]=readest(prm.td,prm.time,'bcp',prm.rcvs{n},prm.dirs.est,...
        prm.fb(4:end),prm.tunit);
    if ~isempty(t)
        data.sats=p.sats; ns=length(p.sats);
        data.bcp{n}=repmat(nan,length(prm.time),ns);
        data.sig{n}=repmat(nan,length(prm.time),ns);
        [tt,i,j]=intersect(prm.time,t);
        for a=arc'
            k=find(a(1)<=tt&tt<=a(2));
            data.bcp{n}(i(k),a(3))=xs(j(k),a(3));
            data.sig{n}(i(k),a(3))=sqrt(vs(j(k),a(3)));
        end
    end
end
set(gcf,'name',['Phase Bias Parameters : ',prm.dirs.est])
set(h,'pointer','arrow');
delete(gut('geth','rcvs'));
gut('newmenu','rcvs','&Receivers',prm.rcvs,prm.rcvs,[mfilename,' cbRcv']);
gut('chkmenu','rcvs',prm.rcv);
data.prm=prm; set(gcf,'userdata',data);

% plot estimated ---------------------------------------------------------------
function plotest(t,x,sig,nsig,ptype,psize)
if ptype==1|ptype==3, plot(t/3600,x,'.','markersize',psize(2)); end
if ptype==2|ptype==3, plot(t/3600,x,'-','linewidth',psize(1)); end
if nsig>0, plot(t/3600,x-nsig*sig,'r:'), plot(t/3600,x+nsig*sig,'r:'), end

% date/time string -------------------------------------------------------------
function str=tstr(td,ts)
str=sprintf('%04d/%02d/%02d %02d:%02d:%02.0f',mjdtocal(td,ts));
