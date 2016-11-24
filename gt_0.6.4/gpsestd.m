function [status,msg]=gpsestd(prm)
%-------------------------------------------------------------------------------
% [system] : GpsTools
% [module] : execute obervation data editor/parameter estimator
% [func]   : execute obervation data editor/parameter estimator
% [argin]  : (prm) = processing parameters struct (see prm_gpsest_def.m)
%                    (default: read from prm_gpsest.m)
% [argout] : status = processing status (1:aborted,0:completed,-1:error)
%            msg    = error message
% [note]   : output files
%            obsc_{rcv}_YYYYDDMMHH.mat : clean observation data
%              epoch = start epoch [year,month,day,hour,min,sec]
%              time  = time vector relative to epoch (sec)
%              sats  = satellite list (cell array)
%              rcv   = station name
%              data  = observation data
%                      data(n,1) : time(n) L1 carrier phase (cycle)
%                      data(n,2) : time(n) L2 carrier phase (cycle)
%                      data(n,3) : time(n) L1 pseudorange (m)
%                      data(n,4) : time(n) L2 pseudorange (m)
%                      data(n,5) : time(n) LC smoothed code (m)
%                      data(n,6) : time(n) L1 smoothed code (m)
%                      data(n,7) : time(n) L2 smoothed code (m)
%              index = observation data index [s,0,f;...]
%                      (s=satellite index,f=1:arc-start,2:arc-end)
%              rstat = receiver states
%                      rstat(n,1)   : time relative to epoch (sec)
%                      rstat(n,2:7) : receiver position/velocity (m,m/sec)
%                      rstat(n,8)   : receiver clock bias (m)
%              azel  = satellite azimath/elevation angle (rad)
%                      azel(n,1) : satellite azimath angle (rad)
%                      azel(n,2) : satellite elevation angle (rad)
%              slip  = cycle-slip positions [t,s,f,...]
%                      (t=time,s=satellite index,f=1:mw,2:gf,3:if,4:td)
%              rpos  = receiver approx. position (m) (ecef)
%              adel  = receiver antenna delta(up/east/north) (m) 
%              atype = receiver antenna type
%              rtype = receiver model type
%            
%            {type}{fb}_{sat|rcv}_YYYYMMDDHH.mat : estimation results
%              epoch = start epoch [year,month,day,hour,min,sec])
%              time  = time vector relative to epoch (sec)
%              sat   = satellite name
%              rcv   = receiver name
%              data  = estimation results
%                      data(n,:) : time(n) estimation results
%              covs  = estimation results variance
%                      covs(n,:) : time(n) estimation results variance
%              prm   = processing parameters struct
%            
%            res{fb}_{rcv}_YYYYMMDDHH.mat : residuals
%              epoch = start epoch [year,month,day,hour,min,sec]
%              time  = time vector relative to epoch (sec)
%              residual = residuals [prefit,postfit,phase bias;...]
%              index = satellite/receiver indexes [s1,s2,r1,r2;...]
%              outl  = outlier flags
%              prm   = processing parameters struct
%            
%            {type}=eph:satellite orbit, clk:clock, zpd:tropspheric parameters,
%                   pos:receiver position, bcp:phase-bias
%            {sat}=satellite name, {rcv}=receiver name, {fb}=f:forward,b:backward
%
% [version]: $Revision: 16 $ $Date: 2008-12-12 15:49:30 +0900 (金, 12 12 2008) $
%            Copyright(c) 2004-2006 by T.Takasu, all rights reserved
% [history]: 04/06/05  0.1  new
%            06/03/23  0.14 restructured
%            07/01/06  0.15 fix bug of saving residual to improper dirs (gt_0.6.3p4)
%                           fix bug of result file name in case of tunit > 24
%            07/05/18  0.16 fix bug of saving debug info (gt_0.6.4)
%            08/11/24  0.17 support session length<=31days (gt_0.6.4)
%                           support time interval>0.01s
%                           support detailed error infomation
%-------------------------------------------------------------------------------
if nargin<1, prm=loadprm('prm_gpsest','prm_gpsest_def'); end
ver=gpstools('version');

status=-1; msg='';
if isempty(prm.sats)|isempty(prm.rcvs), msg='no satellite or receiver'; return, end
if prm.tunit<1|24*366<prm.tunit, msg='estimation unit time out of range'; return, end
if prm.tint<0.01|86400<prm.tint, msg='estimation interval out of range'; return, end
if prm.tover<0|prm.tunit<prm.tover, msg='estimation overlap time out of range'; return, end
if isempty(prm.dirs.est), msg='output directory empty'; return, end

tu=prm.tunit*3600; to=prm.tover*3600; 
[td,ts]=caltomjd(prm.tstart); [tn,te]=caltomjd(prm.tend); te=te+(tn-td)*86400;
if ts<te, prm.ts=ts-to; prm.te=te+to; else msg='no estimation period'; return, end
status=0;

% execute estimation in each processing unit time
for t=(floor(ts/tu):floor(te/tu))*tu
    s.ver=ver; s.ts=now; gt_log('start');
    
    time=max(t,ts)-to:prm.tint:min(t+tu-prm.tint,te)+to;
    gmsg((time(1)-to-prm.ts)/(prm.te-prm.ts));
    if ~prm.dbout
        try
            [status,msg,s]=execest(td,time,prm,ver,s);
        catch
            status=-1;
            gt_log('aborted by error :');
            msg=errlog;
        end
    else
        dbstop if error;
        [status,msg,s]=execest(td,time,prm,ver,s);
        dbclear if error;
    end
    % write processing log
    s.log=gt_log('stop'); s.status=status; s.te=now; gt_writelog(td,time,s,prm);
    if status~=0; break; end
end
if status==0, gmsg(1); end

% execute parameter estimator --------------------------------------------------
function [status,msg,s]=execest(td,time,prm,ver,s)
status=0; msg='';

% initialize common/satelilte parameters
[td,time,prm]=gt_initprm(td,time,prm,1);

% execute observation data editor
if prm.obsedit, [status,msg,s.edit]=gt_obsedit(td,time,prm); end

if prm.gpsest&~status
    % read a priori satellite states
    [state,status,msg,s.ecls]=gt_readstate(td,time,[],prm,1);
    
    if all([prm.satest{:,:},prm.est.erp,prm.est.eco]~=1), m=10; % ppp
    else m=length(prm.rcvs); end
    
    % execute parameter estimator
    for n=1:(length(prm.rcvs)-1)/m+1
        if status, break; else i=n*m-m+1:min(n*m,length(prm.rcvs)); end
        [status,msg,s.est{n}]=estprm(td,time,prm.rcvs(i),state,n,prm);
    end
end

% execute parameter estimator --------------------------------------------------
function [abort,msg,s]=estprm(td,time,rcvs,state,n,prm)
msg=''; s=[]; prm.rcvs=rcvs;

% initialize station parameters
[td,time,prm]=gt_initprm(td,time,prm,2); 

% read a priori station states
[state,abort,msg]=gt_readstate(td,time,state,prm,2);
if abort, return; end

% read clean observation data
[z,tz,iz,arc,antp,state,abort,s.obss]=gt_readobs(td,time,state,prm);
if abort|isempty(z), return; end

% initialize state variables
[x,ix,P,sig,prn]=gt_initstate(prm);

% pass-1
[x,P,z,arc,abort,s.est1]=...
    estpass(td,time,x,ix,P,z,tz,iz,antp,state,arc,sig,prn,[n,1],1,prm);

if prm.backward&~abort % pass-2
    [P,arc]=resetcov(ix,P,arc,sig,prm);
    [x,P,z,arc,abort,s.est2]=...
        estpass(td,time,x,ix,P,z,tz,iz,antp,state,arc,sig,prn,[n,2],2,prm);
end
if prm.iteration&~abort % pass-3
    [P,arc]=resetcov(ix,P,arc,sig,prm);
    [x,P,z,arc,abort,s.est3]=...
        estpass(td,time,x,ix,P,z,tz,iz,antp,state,arc,sig,prn,[n,3],1,prm);
end

% estimation pass --------------------------------------------------------------
function [x,P,z,arc,abort,stats]=estpass(td,time,x,ix,P,z,tz,iz,antp,state,arc,...
                                         sig,prn,pass,fb,prm)
x0=x; P0=P; arc0=arc; abort=0; stats={}; nz=size(z,1); if nz==0, return, end
[tn,ts]=caltomjd(prm.tstart); to=(td-tn)*86400;
[t,i,j]=intersect(time,tz(:,2)); j=[1,j+1];

for pasn=1:10
    xs=repmat(nan,length(time),length(x)); ps=xs;
    phs=zeros(length(prm.sats),length(prm.rcvs));
    res=zeros(nz,3); if prm.dbout, dbg=zeros(nz,17); else dbg=[]; end
    
    if fb==1, nn=1:length(t); p=1; else nn=length(t):-1:1; p=-1; end
    for n=nn
        q=(to+t(n)-prm.ts)/(prm.te-prm.ts);
        if ticks(2)&gmsg(q), abort=1; break; end
        
        % measurement update of states
        k=j(n):j(n+1)-1;
        [x,P,arc,phs,iz(k,3),res(k,:),iv,ds]=gt_udmeas(td,t(n),x,ix,P,arc,...
            phs,z(k,:),tz(k,1),iz(k,:),antp,state(i(n)),sig,[pass,pasn],fb,prm);
        
        xs(i(n),:)=x'; ps(i(n),:)=diag(P)'; x(iv)=nan; P(iv,iv)=nan;
        if prm.dbout, dbg(k,:)=ds; end
        
        if n~=nn(end) % temporal update of states
            [x,P]=gt_udstate(td,t(n),t(n+p),x,ix,P,state(i(n)),prn,prm);
        end
    end
    % quality control by statistics
    gmsg(q); gmsg('pass-%d-%d-%d : quality checking',pass,pasn);
    [xs,xz,xa,stats]=gt_qc(td,time,xs,ps,ix,tz(:,2),iz,res,arc,prm,pass,abort);
    z(xz,:)=nan; arc(xa,1:2)=nan;
    
    if abort|~prm.reests|isempty(xz), break; end
    
    x=x0; P=P0; arc0(xa,1:2)=nan; arc=arc0; % rollback
end
% quality control/alignment of estimated clock
if ~abort, xs=gt_alignclk(td,time,xs,ix,stats,prm); end

% save estimation results/residuals
if ~abort
    saveest(td,time,ix,xs,ps,arc,fb,prm,to==0);
    saveres(td,time(1),tz(:,2),iz,res,fb,prm,to==0);
end
if ~isempty(dbg), savedbg(td,tz(:,2),iz,dbg,fb,prm); end

% reset covariences ------------------------------------------------------------
function [P,arc]=resetcov(ix,P,arc,sig,prm)
if prm.sdfact>0
    P=P*prm.sdfact.^2; if ~isempty(arc), arc(:,8)=arc(:,8)*prm.sdfact; end
else
    P=diag(sig.^2); if ~isempty(arc), arc(:,8)=prm.sig.arcn(arc(:,4),1); end
end

% save estimation results ------------------------------------------------------
function saveest(td,time,ix,xs,ps,arc,fb,prm,opt)
t=time(1)+prm.tover*3600; if opt, tu=prm.tunit*3600; t=floor(t/tu)*tu; end
epoch=mjdtocal(td,t); time=time(:)-t;
arc(:,1:2)=arc(:,1:2)-t;
for n=1:length(prm.sats)
    if ~isempty(arc), arcs=arc(arc(:,3)==n,[1:2,4]); else arcs=[]; end
    i=[ix.sato{n},ix.sats{n}]; j=ix.satc{n};
    saveestf('eph',fb,epoch,time,prm.sats{n},'',xs(:,i),ps(:,i),arcs,prm);
    saveestf('clk',fb,epoch,time,prm.sats{n},'',xs(:,j),ps(:,j),arcs,prm);
end
for n=1:length(prm.rcvs)
    if ~isempty(arc), arcr=arc(arc(:,4)==n,[1:2,3]); else arcr=[]; end
    i=ix.rcvc{n}; j=[ix.rcvz{n},ix.rcvg{n}]; k=ix.rcvp{n}; l=[ix.arcn{:,n}]';
    saveestf('clk',fb,epoch,time,'',prm.rcvs{n},xs(:,i),ps(:,i),arcr,prm);
    saveestf('zpd',fb,epoch,time,'',prm.rcvs{n},xs(:,j),ps(:,j),arcr,prm);
    saveestf('pos',fb,epoch,time,'',prm.rcvs{n},xs(:,k),ps(:,k),arcr,prm);
    saveestf('bcp',fb,epoch,time,'',prm.rcvs{n},xs(:,l),ps(:,l),arcr,prm);
end
i=ix.erp;
if ~isempty(i)
    f=repmat([3600*180/pi,3600*180/pi,1],size(xs,1),1);
    saveestf('erp',fb,epoch,time,'','',xs(:,i).*f,ps(:,i).*f.^2,arc,prm);
end
i=ix.eco; saveestf('eco',fb,epoch,time,'','',xs(:,i),ps(:,i),arc,prm);

% save estimation result to file -----------------------------------------------
function saveestf(type,fb,epoch,time,sat,rcv,data,covs,arc,prm)
if ~prm.estout|isempty(data), return, end
s='fb'; if isempty(sat), name=rcv; else name=sat; end
file=sprintf('%s%s_%s_%04d%02d%02d%02d.mat',type,s(fb),name,epoch(1:4));
file=gfilepath(prm.dirs.est,file,epoch,name,1);
gmsg('saving : %s',file);
save(file,'epoch','time','sat','rcv','data','covs','arc','prm')

% save resuduals ---------------------------------------------------------------
function saveres(td,ts,tz,iz,res,fb,prm,opt)
if ~prm.statout|isempty(tz), return, end
t=ts+prm.tover*3600; if opt, tu=prm.tunit*3600; t=floor(t/tu)*tu; end
epoch=mjdtocal(td,t); s='fb';
for n=1:length(prm.rcvs)
    i=find(iz(:,2)==n); time=tz(i)-t; index=iz(i,[1,1,2,2]); index(:,[2,4])=0;
    residual=res(i,:); outl=iz(i,3)==2;
    file=sprintf('res%s_%s_%04d%02d%02d%02d.mat',s(fb),prm.rcvs{n},epoch(1:4));
    file=gfilepath(prm.dirs.est,file,epoch,prm.rcvs{n},1);
    gmsg('saving : %s',file);
    save(file,'epoch','time','index','residual','outl','prm');
end

% save debug information -------------------------------------------------------
function savedbg(td,tz,iz,dbg,fb,prm)
s='fb'; [time,i]=sort(tz); index=iz(i,:); dbg=dbg(i,:); epoch=datevec(now);
if ~isempty(prm.rcvs),name=prm.rcvs{1}; else name=''; end
file=sprintf('dbg%s_%04d%02d%02d%02d%02d%02.0f.mat',s(fb),epoch);
file=gfilepath(prm.dirs.est,file,prm.tstart,name,1);
gmsg('saving : %s',file);
save(file,'td','time','index','dbg','prm');

% program error log ------------------------------------------------------------
function msg=errlog
err=lasterror;
gt_log('  identifier: %s',err.identifier);
gt_log('  message   : %s',err.message);
for n=1:length(err.stack)
    gt_log('  stack info: (%d) L%d %s %s',n,err.stack(n).line,...
           err.stack(n).name,err.stack(n).file);
end
gt_log('  workspace :');
if ~isempty(err.stack)
    msg=sprintf('%s at L%d in %s',err.identifier,err.stack(1).line,...
                err.stack(1).file);
else
    msg=err.message;
end

% check tick -------------------------------------------------------------------
function f=ticks(n)
persistent tp, ts=now; if 86400*(ts-tp)<1/n, f=0; else f=1; tp=ts; end
