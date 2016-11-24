function [td,time,prm]=gt_initprm(td,time,prm,opt)
%-------------------------------------------------------------------------------
% [system] : GpsTools
% [module] : initialize processing parameters
% [func]   : initialize processing parameters
% [argin]  : td,time = input date(mjd-gpst),start time (sec)
%            prm     = input processing parameters (see prm_gpsest_def.m)
%            (opt)   = options (1:only satellite,2:only station)
%                    (default:0)
% [argout] : td,time = output date(mjd-gpst),start time (sec)
%            prm     = output processing parameters
% [note]   :
% [version]: $Revision: 12 $ $Date: 06/07/21 1:24 $
%            Copyright(c) 2004-2006 by T.Takasu, all rights reserved
% [history]: 06/02/14  0.1  separated from gpsestd.m
%            06/06/25  0.2  add error message
%-------------------------------------------------------------------------------
if nargin<4, opt=0; end
if opt~=2, [td,time,prm]=inits(td,time,prm); end
if opt~=1, prm=initr(td,time,prm); end

% initialize common/satellite parameters ---------------------------------------
function [td,time,prm]=inits(td,time,prm)

gmsg('initializing common/satellite parameters');

dd=floor(time(1)/86400);
td=td+dd; time=time-86400*dd;
tstart=mjdtocal(td,time(1));

prm.f1=prm.f1*1E9;
prm.f2=prm.f2*1E9;
prm.cif=[prm.f1^2;-prm.f2^2]/(prm.f1^2-prm.f2^2);
prm.elmin=prm.elmin*pi/180;
prm.elmax=prm.elmax*pi/180;
prm.obs.elmin=prm.obs.elmin*pi/180;

% extend file path
prm.dirs.est=gfilepath('',prm.dirs.est,tstart);
if ~isempty(prm.dirs.obc), prm.dirs.obc =gfilepath('',prm.dirs.obc,tstart); else prm.dirs.obc=prm.dirs.est; end
if ~isempty(prm.dirs.obs), prm.dirs.obs =gfilepath('',prm.dirs.obs, []); else prm.dirs.obs =prm.dirs.est; end
if ~isempty(prm.dirs.nav), prm.dirs.nav =gfilepath('',prm.dirs.nav, []); else prm.dirs.nav =prm.dirs.est; end
if ~isempty(prm.dirs.eph), prm.dirs.eph =gfilepath('',prm.dirs.eph, []); else prm.dirs.eph =prm.dirs.est; end
if ~isempty(prm.dirs.clk), prm.dirs.clk =gfilepath('',prm.dirs.clk, []); else prm.dirs.clk =prm.dirs.est; end
if ~isempty(prm.dirs.pos), prm.dirs.pos =gfilepath('',prm.dirs.pos, []); else prm.dirs.pos =prm.dirs.est; end
if ~isempty(prm.dirs.erp), prm.dirs.erp =gfilepath('',prm.dirs.erp, []); else prm.dirs.erp =prm.dirs.est; end
if ~isempty(prm.dirs.trop),prm.dirs.trop=gfilepath('',prm.dirs.trop,[]); else prm.dirs.trop=prm.dirs.est; end
if ~isempty(prm.dirs.ion), prm.dirs.ion =gfilepath('',prm.dirs.ion, []); else prm.dirs.ion =prm.dirs.est; end
if ~isempty(prm.dirs.inp), prm.dirs.inp =gfilepath('',prm.dirs.inp, []); else prm.dirs.inp =prm.dirs.est; end
if ~isempty(prm.dirs.inpr),prm.dirs.inpr=gfilepath('',prm.dirs.inpr,[]); else prm.dirs.inpr=prm.dirs.est; end
prm.satsrpf =gfilepath('',prm.satsrpf, tstart);
prm.rcvposf =gfilepath('',prm.rcvposf, tstart);
prm.pcv     =gfilepath('',prm.pcv,     tstart);
prm.satpcv  =gfilepath('',prm.satpcv,  tstart);
prm.oload   =gfilepath('',prm.oload,   tstart);
prm.mpc     =gfilepath('',prm.mpc,     tstart);
prm.mpp     =gfilepath('',prm.mpp,     tstart);
prm.exsatrcv=gfilepath('',prm.exsatrcv,tstart);

% read satellite parameters
prm.satprm=prm_gpssats;

% set eclipse parameters
prm.satblock=zeros(length(prm.sats),1);
prm.eclnsigs=repmat(prm.eclnsig(1),length(prm.sats),1);
prm.eclnprns=repmat(prm.eclnprn(1),length(prm.sats),1);
if ~isempty(prm.satprm)
    [s,i,j]=intersect(prm.sats,prm.satprm(:,1));
    i=i([prm.satprm{j,2}]==1); % Block IIR
    prm.satblock(i)=1;
    prm.eclnsigs(i)=prm.eclnsig(2);
    prm.eclnprns(i)=prm.eclnprn(2);
end
% read satellite orbit parameters
prm.srpprms={};
[dirs,file,ext]=fileparts(prm.satsrpf);
if exist(prm.satsrpf)
    wd=pwd; cd(dirs); prm.srpprms=feval(file); cd(wd);
elseif ~isempty(prm.satsrpf)
    gt_log('no sat orbit parameter  : %s',prm.satsrpf);
end
% read satellite pcv parameters
prm.spcvprms={};
[dirs,file,pcvext]=fileparts(prm.satpcv);
if exist(prm.satpcv)
    if strcmp(pcvext,'.m')
        wd=pwd; cd(dirs); prm.spcvprms=feval(file); cd(wd);
    elseif ~strcmp(pcvext,'.atx')
        gt_log('no valid sat ant pcv    : %s',prm.satpcv);
    end
elseif ~isempty(prm.satpcv)
    gt_log('no sat antenna pcv file : %s',prm.satpcv);
end
prm.ants=zeros(3+16,length(prm.sats));
for n=1:length(prm.sats)
    if isempty(prm.satest), error('@no sat parameter settings'), end
    i=find(strcmp(prm.sats{n},prm.satest(:,1)));
    if isempty(i), i=find(strcmp('ALL',prm.satest(:,1))); end
    if isempty(i), error(['@no sat parameter setting : ',prm.sats{n}]), end
    prm.est.sato(n,1)=prm.satest{i,2};
    prm.est.sats(n,1)=prm.satest{i,3};
    prm.est.satc(n,1)=prm.satest{i,4};
    
    if isempty(prm.satsig), error('@no sat a priori std deviation settings'), end
    i=find(strcmp(prm.sats{n},prm.satsig(:,1)));
    if isempty(i), i=find(strcmp('ALL',prm.satsig(:,1))); end
    if isempty(i), error(['@no sat a priori std deviation setting : ',prm.sats{n}]), end
    prm.sig.sato(n,:)=[prm.satsig{i,[2,2,2,3,3,3]}];
    prm.sig.sats(n,:)=[prm.satsig{i,[4,5,6,7,7,7]}];
    prm.sig.satc(n,:)=[prm.satsig{i,[8,9]}];
    
    if isempty(prm.satprn), error('@no sat process noise settings'), end
    i=find(strcmp(prm.sats{n},prm.satprn(:,1)));
    if isempty(i), i=find(strcmp('ALL',prm.satprn(:,1))); end
    if isempty(i), error(['@no sat process noise setting : ',prm.sats{n}]), end
    prm.prn.sato(n,:)=[prm.satprn{i,[2,2,2,3,3,3]}];
    prm.prn.sats(n,:)=[prm.satprn{i,[4,5,6,7,7,7]}];
    prm.prn.satc(n,:)=[prm.satprn{i,[8,9]}];
    
    if isempty(prm.satprm), error('@no sat parameters'), end
    i=find(strcmp(prm.sats{n},prm.satprm(:,1)));
    if isempty(i), i=find(strcmp('ALL',prm.satprm(:,1))); end
    if isempty(i), error(['@no sat parameter : ',prm.sats{n}]), end
    prm.sat.srp(n,:)=[prm.satprm{i,2:end}];
    
    % satellite antenna offset
    if strcmp(pcvext,'.atx') % antex pcv file
        [apc1,apc2,apv1,apv2]=readpcv('',prm.sats{n},prm.satpcv,td+time(1)/86400);
        prm.ants(:,n)=[apc1,apv1(:,1:16)]';
    else
        i=find(strcmp(prm.sats{n},prm.satprm(:,1)));
        if ~isempty(i)
            if prm.satprm{i,2}==0, prm.ants(1:3,n)=[0.279;0.000;1.023]; end
        else
            gt_log('no satellite parameters     : %s',prm.sats{n});
        end
        % satellite antenna pcv parameter
        if ~isempty(prm.spcvprms)
            i=find(strcmp(prm.sats{n},prm.spcvprms(:,1)));
            if ~isempty(i)
                prm.ants(4:end,n)=-[prm.spcvprms{i,2:end}]';
            else
                gt_log('no satellite antenna pcv      : %s file=%s',prm.sats{n},prm.satpcv);
            end
        end
    end
    % satellite srp parameter
    if ~isempty(prm.srpprms)
        i=find(strcmp(prm.sats{n},prm.srpprms(:,1)));
        if ~isempty(i)
            prm.sig.sats(n,:)=[prm.srpprms{i,8:13}];
            prm.prn.sats(n,:)=[prm.srpprms{i,14:19}]/10;
        else
            gt_log('no satellite srp params : %s file=%s',prm.sats{n},prm.satsrpf);
        end
    end
end
% read initial earth rotation parameter
if prm.est.erp>=1&~isempty(prm.initerp)
    prm.erp=readerp(tutc(td,time(1)),prm.dirs.erp,prm.initerp);
else
    prm.erp=[0,0,0,0,0];
end
prm.sig.erp=prm.sig.erp([1,1,2]);
prm.sig.eco=prm.sig.eco([1,1,1,2]);
prm.prn.erp=prm.prn.erp([1,1,2]);
prm.prn.eco=prm.prn.eco([1,1,1,2]);

% read excluded satellite/receiver
if ~isempty(prm.exsatrcv)
    [dirs,file,ext]=fileparts(prm.exsatrcv);
    if exist(prm.exsatrcv)
        wd=pwd; cd(dirs); exclude=feval(file); cd(wd);
        prm.exclude=[prm.exclude;exclude];
    elseif ~isempty(prm.exsatrcv)
        gt_log('no excluded sat/rcv file: %s',prm.exsatrcv);
    end
end

% initialize receiver parameters -----------------------------------------------
function prm=initr(td,time,prm)

gmsg('initializing receiver parameters');

if isempty(prm.tpos), prm.tp=[]; else [prm.tp,t]=caltomjd(prm.tpos); end
prm.mpcc=[]; prm.mpcs=[];

% baselines
prm.ircv=[];
for n=1:size(prm.baseline,1)
    i=find(strcmp(prm.baseline{n,1},prm.rcvs));
    if ~isempty(i)
        j=find(strcmp(prm.baseline{n,2},prm.rcvs));
        if ~isempty(j), prm.ircv=[prm.ircv;i,j];
        elseif strcmp(prm.baseline{n,2},'ALL')
            for j=1:length(prm.rcvs), if i~=j, prm.ircv=[prm.ircv;i,j]; end, end
        end
    end
end
% read receiver positions
prm.rcvposs={};
[dirs,file,ext]=fileparts(prm.rcvposf);
if exist(prm.rcvposf)
    wd=pwd; cd(dirs); prm.rcvposs=feval(file); cd(wd);
elseif ~isempty(prm.rcvposf)
    gt_log('no receiver pos file    : %s',prm.rcvposf);
end
for n=1:length(prm.rcvs)
    if isempty(prm.rcvest), error('@no rcv parameter settings'), end
    i=find(strcmp(prm.rcvs{n},prm.rcvest(:,1)));
    if isempty(i), i=find(strcmp('ALL',prm.rcvest(:,1))); end
    if isempty(i), error(['@no rcv parameter setting : ',prm.rcvs{n}]), end
    prm.est.rcvc(n,1)=prm.rcvest{i,2};
    prm.est.rcvz(n,1)=prm.rcvest{i,3};
    prm.est.rcvg(n,1)=prm.rcvest{i,3}&prm.trgmodel>=1;
    prm.est.rcvp(n,1)=prm.rcvest{i,4};
    prm.est.rcvb(n,1)=prm.rcvest{i,5};
    if strcmp(prm.rcvs{n},prm.clkref), prm.est.rcvc(n,1)=0; end % reference clock
    
    if isempty(prm.obssig), error('@no measurement noise settings'), end
    i=find(strcmp(prm.rcvs{n},prm.obssig(:,1)));
    if isempty(i), i=find(strcmp('ALL',prm.obssig(:,1))); end
    if isempty(i), error(['@no measurement noise setting : ',prm.rcvs{n}]), end
    prm.sig.obs(n,1)=prm.obssig{i,2};
    
    if isempty(prm.rcvsig), error('@no rcv a priori std deviation settings'), end
    i=find(strcmp(prm.rcvs{n},prm.rcvsig(:,1)));
    if isempty(i), i=find(strcmp('ALL',prm.rcvsig(:,1))); end
    if isempty(i), error(['@no rcv a priori std deviation setting : ',prm.rcvs{n}]), end
    prm.sig.rcvc(n,:)=[prm.rcvsig{i,2:3}];
    prm.sig.rcvz(n,:)=[prm.rcvsig{i,4:5}];
    prm.sig.rcvg(n,:)=[prm.rcvsig{i,[6,6,6,6,6]}];
    prm.sig.rcvp(n,1:3)=[prm.rcvsig{i,7:9}];
    if prm.rposmodel>=2, prm.sig.rcvp(n,4:6)=[prm.rcvsig{i,7:9}]/30; end
    prm.sig.arcn(n,:)=[prm.rcvsig{i,10}];
    
    if isempty(prm.rcvsig), error('@no rcv process noise settings'), end
    i=find(strcmp(prm.rcvs{n},prm.rcvprn(:,1)));
    if isempty(i), i=find(strcmp('ALL',prm.rcvprn(:,1))); end
    if isempty(i), error(['@no rcv process noise setting : ',prm.rcvs{n}]), end
    prm.prn.rcvc(n,:)=[prm.rcvprn{i,2:3}];
    prm.prn.rcvz(n,:)=[prm.rcvprn{i,4:5}];
    prm.prn.rcvg(n,:)=[prm.rcvprn{i,[6,6,6,6,6]}];
    prm.prn.rcvp(n,:)=[prm.rcvprn{i,7:9}];
    prm.prn.arcn(n,:)=[prm.rcvprn{i,10}];
    
    % read phase multipath profile
    if ~isempty(prm.mpp)
        file=gfilepath('',prm.mpp,[],prm.rcvs{n});
        mpcc=[]; mpcs=[];
        if exist(file)
            load(file)
        else
            gt_log('no phase multipath file : %s',file);
        end
        if ~isempty(mpcc)&~isempty(mpcs)
            if isempty(prm.mpcc)
                prm.mpcc=zeros([size(mpcc),length(prm.rcvs)]); prm.mpcs=prm.mpcc;
                prm.mpcc(:,:,n)=mpcc; prm.mpcs(:,:,n)=mpcs;
            elseif all(size(mpcc)==size(prm.mpcc(:,:,1)))&...
                   all(size(mpcs)==size(prm.mpcs(:,:,1)))
                prm.mpcc(:,:,n)=mpcc; prm.mpcs(:,:,n)=mpcs;
            else
                gt_log('multipath imcompatible : %s',file);
            end
        end
    end
end
% read ocean loading parameters
[prm.odisp,prm.ophas]=readoload(prm.rcvs,prm.oload);

% gpst->mjd-utc ----------------------------------------------------------------
function [t,utc_tai]=tutc(td,t)
utc_tai=prm_utc_tai(td+t/86400,1);
t=td+(t+19+utc_tai)/86400;
