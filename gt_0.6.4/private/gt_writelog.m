function gt_writelog(td,time,stats,prm)
%-------------------------------------------------------------------------------
% [system] : GpsTools
% [module] : write processing log file
% [func]   : write processing log file
% [argin]  : td,time = day(mjd),time vector(sec)
%            stats   = processing status/statistics
%            prm     = processing parameters
% [argout] : none
% [note]   :
% [version]: $Revision: 16 $ $Date: 2008-12-12 15:49:30 +0900 (金, 12 12 2008) $
%            Copyright(c) 2004-2008 by T.Takasu, all rights reserved
% [history]: 06/02/14  0.1  separated from gpsestd.m
%            08/12/02  0.2  fix bug on writing clock-jump log
%                           fix bug on writing eclipse log
%                           support new satclk, rcvpos
%-------------------------------------------------------------------------------
dv=datevec(stats.ts); sstr={'error','completed','aborted'};
file=sprintf('gpsestd_%04d%02d%02d%02d%02d%02.0f.log',dv);
file=gfilepath(prm.dirs.est,file,mjdtocal(td),'',1);
f=fopen(file,'wt');
if f==0, disp(['warning : log file open error : ',file]), return, end
fprintf(f,'*****     GpsTools ver.%s : OBS DATA EDITOR/PARAMETER ESTIMATOR LOG     *****\n\n',stats.ver);
fprintf(f,'start/end time                      : %s\n',pstr(datevec(stats.ts),datevec(stats.te)));
fprintf(f,'processing time                     : %s\n',tstr(datevec(stats.te-stats.ts)));
fprintf(f,'processing status                   : %s\n',sstr{stats.status+2});
fprintf(f,'\n');
fprintf(f,'processing parameters :\n\n');
writeprm1(f,td,time,prm);
if prm.obsedit, writeprm2(f,prm); end
if prm.gpsest,  writeprm3(f,prm); end
if isfield(stats,'edit'), writeedit(f,stats.edit,prm); end
obss={}; est1={}; est2={}; est3={};
if isfield(stats,'est')
    for n=1:length(stats.est)
        s=stats.est{n};
        if isfield(s,'obss'), obss={obss{:},s.obss{:}}; end
        if isfield(s,'est1'), est1={est1{:},s.est1{:}}; end
        if isfield(s,'est2'), est2={est2{:},s.est2{:}}; end
        if isfield(s,'est3'), est3={est3{:},s.est3{:}}; end
    end
end
if ~isempty(obss), writeobss(f,obss,prm); end
if isfield(stats,'ecls'), writeecls(f,stats.ecls,prm); end
if ~isempty(est1), writeests(f,est1,prm,1); end
if ~isempty(est2), writeests(f,est2,prm,2); end
if ~isempty(est3), writeests(f,est3,prm,3); end
fprintf(f,'processing logs :\n\n');
fprintl(f,stats.log);
fclose(f);

% common processing parameters -------------------------------------------------
function writeprm1(f,td,time,prm)
if isempty(time), time=0; end
ts=mjdtocal(td,time(1));
te=mjdtocal(td,time(end));
fprintf(f,'estimation start time               : %s GPST (GPSW%d-%d)\n',tstr(ts),gpswd(ts));
fprintf(f,'estimation end time                 : %s GPST (GPSW%d-%d)\n',tstr(te),gpswd(te));
fprintf(f,'estimation interval                 : %.1f sec\n',prm.tint);
fprintf(f,'estimation unit time                : %.1f hr\n',prm.tunit);
fprintf(f,'estimation overlap time             : %.1f hr\n',prm.tover);
fprintf(f,'\n');
fprintf(f,'satellites (%d)\n  ',length(prm.sats));
for n=1:length(prm.sats)
    if mod(n,12)~=0, s=' '; elseif n<length(prm.sats), s='\n  '; end
    fprintf(f,'%-5s',prm.sats{n}); fprintf(f,s);
end
fprintf(f,'\n\n');
fprintf(f,'receivers (%d)\n  ',length(prm.rcvs));
for n=1:length(prm.rcvs)
    if mod(n,10)~=0, s=' '; elseif n<length(prm.rcvs), s='\n  '; end
    fprintf(f,'%-6s ',prm.rcvs{n}); fprintf(f,s);
end
fprintf(f,'\n\n');
fprintf(f,'execute observation data editor     : %s\n',onoff(prm.obsedit));
fprintf(f,'execute parameter estimator         : %s\n',exepass(prm));
fprintf(f,'raw observation data                : %s\n',prm.obs.src);
fprintf(f,'navigation messages                 : %s\n',prm.src.nav);
fprintf(f,'clean observation data directory    : %s\n',prm.dirs.obc);
fprintf(f,'output directory                    : %s\n',prm.dirs.est);
fprintf(f,'data directories\n');
fprintf(f,'  raw observation data              : %s\n',prm.dirs.obs);
fprintf(f,'  navigation messages               : %s\n',prm.dirs.nav);
fprintf(f,'  satellite ephemerides             : %s\n',prm.dirs.eph);
fprintf(f,'  satellite/receiver clocks         : %s\n',prm.dirs.clk);
fprintf(f,'  receiver positions                : %s\n',prm.dirs.pos);
fprintf(f,'  earth rotation parameters         : %s\n',prm.dirs.erp);
fprintf(f,'  tropospheric parameters           : %s\n',prm.dirs.trop);
fprintf(f,'  input estimated sats. parameters  : %s\n',prm.dirs.inp);
fprintf(f,'  input estimated rcvs. parameters  : %s\n',prm.dirs.inpr);
fprintf(f,'data files\n');
fprintf(f,'  satellite antenna pcv parameter   : %s\n',prm.satpcv);
fprintf(f,'  receiver antenna pcv parameter    : %s\n',prm.pcv);
fprintf(f,'  receiver position                 : %s\n',prm.rcvposf);
fprintf(f,'  ocean loading parameter           : %s\n',prm.oload);
fprintf(f,'  excluded satellites/receviers     : %s\n',prm.exsatrcv);
fprintf(f,'  phase multipath profile           : %s\n',prm.mpc);
fprintf(f,'execution environment\n');
fprintf(f,'  matlab version                    : %s\n',version);
fprintf(f,'  operating system                  : %s\n',computer);
fprintf(f,'  java virtual machine              : %s\n',version('-java'));
fprintf(f,'\n');

% observation data editor parameters ---------------------------------------------------
function writeprm2(f,prm)
fprintf(f,'observation data editor settings\n');
fprintf(f,'  L1 frequency                      : %.6f GHz\n',prm.f1);
fprintf(f,'  L2 frequency                      : %.6f GHz\n',prm.f2);
fprintf(f,'  observation time tolerance        : %.3f sec\n',prm.ttol);
fprintf(f,'  exclude unhealthy satellites      : %s\n',onoff(prm.exuhsat));
fprintf(f,'  min elevation angle               : %.1f deg\n',prm.obs.elmin);
fprintf(f,'  max arc gap                       : %.1f sec\n',prm.obs.gapmax);
fprintf(f,'  min arc length                    : %.1f sec\n',prm.obs.arcmin);
fprintf(f,'  min arc points                    : %.1f points\n',prm.obs.pntmin);
fprintf(f,'  noise of Melbourne-Wubbena code   : %.2f m\n',prm.obs.sigmw);
fprintf(f,'  noise of ion-free code            : %.2f m\n',prm.obs.sigif);
fprintf(f,'  cycle-slip threshold of LG        : %.2f m\n',prm.obs.dgmax);
fprintf(f,'  cycle-slip threshold of WL        : %.2f cycle\n',prm.obs.dwl);
fprintf(f,'  cycle-slip threshold of MP1       : %.2f m\n',prm.obs.dmp1);
fprintf(f,'  cycle-slip threshold of MP2       : %.2f m\n',prm.obs.dmp2);
fprintf(f,'  outlier threshold of WL           : %.1f sigmas\n',prm.obs.outl);
fprintf(f,'  ionos smoothing window size       : %.0f points\n',prm.ionwind);
fprintf(f,'  moving average window size (1)    : %.0f points\n',prm.obs.wind(1));
fprintf(f,'  moving average window size (2)    : %.0f points\n',prm.obs.wind(2));
fprintf(f,'  slip thres. elevation weighting   : %s\n',onoff(prm.obs.elwe));
fprintf(f,'  repair steered clock jump         : %s\n',obsopt1(prm.obs.clkrep));
fprintf(f,'  repair cycle-slip                 : %s\n',onoff(prm.obs.reps));
fprintf(f,'  no. of points for poly. fitting   : %.0f\n',prm.obs.npnt);
fprintf(f,'  degree of polynomial fitting      : %.0f\n',prm.obs.nmax);
fprintf(f,'  code multipath profile            : %s\n',prm.mpc);
fprintf(f,'  separate arc at eclipse in/out    : %s\n',onoff(prm.obs.separc));
fprintf(f,'\n');

% parameter estimator parameters -----------------------------------------------
function writeprm3(f,prm)
fprintf(f,'parameter estimator settings\n');
fprintf(f,'  estimation strategy               : %s\n',prm.omodel);
fprintf(f,'  satellite clock model             : %s\n',modelclk(prm.sclkmodel));
fprintf(f,'  receiver clock model              : %s\n',modelclk(prm.rclkmodel));
fprintf(f,'  receiver position model           : %s\n',modelpos(prm.rposmodel));
fprintf(f,'  clock reference                   : %s\n',prm.clkref);
fprintf(f,'estimation/measurement model settings\n');
fprintf(f,'  min elevation angle               : %.1f deg\n',prm.elmin);
fprintf(f,'  max elevation angle               : %.1f deg\n',prm.elmax);
fprintf(f,'  tropospheric a priori model       : %s\n',prm.trop);
fprintf(f,'  tropospheric mapping function     : %s\n',prm.mapf);
fprintf(f,'  tropospheric stochastic model     : %s\n',modelzpd(prm.zpdmodel));
fprintf(f,'  tropospheric gradient model       : %s\n',modeltrg(prm.trgmodel));
fprintf(f,'  meas. noise elevation angle factor: %.1f\n',prm.sigweight);
fprintf(f,'  meas. noise eclipse factor        : %4.1f %4.1f\n',prm.eclnsig);
fprintf(f,'  process noise factor in post-ecl. : %4.1f %4.1f\n',prm.eclnprn);
fprintf(f,'  max post-eclipse maneuver time    : %.1f sec\n',prm.ecltime);
fprintf(f,'  reset time of observation outage  : %.1f sec\n',prm.ogapmax);
fprintf(f,'  std. dev. factor of initial states: %.1f\n',prm.sdfact);
fprintf(f,'  std. dev. factor of pre-det states: %.1f\n',prm.sdfact2);
fprintf(f,'  max count of filter iteration     : %.0f\n',prm.maxiter);
fprintf(f,'  min count of obs for estimation   : %.0f\n',prm.minobs);
fprintf(f,'  nutation model                    : %s\n',prm.nutmodel);
fprintf(f,'  filter                            : %s\n',prm.filter);
fprintf(f,'  clock correction by phase-bias    : %s\n',onoff(prm.clkcorr));
fprintf(f,'  qc clock wrt reference            : %s\n',prm.clkqc);
fprintf(f,'  align clock to reference          : %s\n',prm.clkalign);
fprintf(f,'  receiver attitude model           : %s\n',modelatt(prm.rattmodel));
fprintf(f,'  measurement corrections\n');
fprintf(f,'    satellite antenna offsets       : %s\n',onoff(prm.corrf(1)));
fprintf(f,'    receiver antenna offsets        : %s\n',onoff(prm.corrf(2)));
fprintf(f,'    relativistic effects            : %s\n',coropt1(prm.corrf(3)));
fprintf(f,'    phase windup                    : %s\n',onoff(prm.corrf(4)));
fprintf(f,'    phase multipath                 : %s\n',onoff(prm.corrf(5)));
fprintf(f,'  station displacement corrections\n');
fprintf(f,'    solid earth tides               : %s\n',onoff(prm.sitedisp(1)));
fprintf(f,'    ocean loading                   : %s\n',onoff(prm.sitedisp(2)));
fprintf(f,'    pole tides                      : %s\n',onoff(prm.sitedisp(3)));
fprintf(f,'    eliminate permanent deformation : %s\n',onoff(prm.sitedisp(4)));
fprintf(f,'  erp variation correction          : %s\n',onoff(prm.erpvar));
fprintf(f,'  meteorological parameters         : %s\n',prm.metsrc);
fprintf(f,'  default mean sea level pressure   : %.1f hPa\n',prm.metprm(1));
fprintf(f,'  default temperature/humidity      : %.1f C %.1f %%\n',prm.metprm(2:3));
if any([prm.satest{:,2}]==1)
    fprintf(f,'satellite orbit model settings\n');
    fprintf(f,'  max degree of geopotential        : %.0f\n',prm.obt.g_nmax);
    fprintf(f,'  solar/planetary potentials        : %s\n',prm.obt.p_plgrv);
    fprintf(f,'  solar radiation pressure model    : %s\n',prm.obt.p_solarpr);
    fprintf(f,'  eclipse model                     : %s\n',prm.obt.p_eclipse);
    fprintf(f,'  atmospheric drag model            : %s\n','');
    fprintf(f,'  earth tides                       : %s\n',prm.obt.p_tidal);
    fprintf(f,'  relativistic effects              : %s\n',prm.obt.p_relativ);
    fprintf(f,'  thrust forces                     : %s\n',prm.obt.p_deltav);
    fprintf(f,'  integration step                  : %.1f sec\n',prm.obt.tstep);
end
fprintf(f,'quality control settings\n');
fprintf(f,'  re-estimate without exclusions    : %s\n',onoff(prm.reests));
fprintf(f,'  outlier threshold of prefit res.  : %.1f sigmas\n',prm.outlp);
fprintf(f,'  outlier threshold of postfit res. : %.1f sigmas\n',prm.outlf);
fprintf(f,'  arc max outlier rate              : %.0f %%\n',prm.varc.pout*100);
fprintf(f,'  arc max postfit residual rms      : %.2f m\n',prm.varc.rmsf);
fprintf(f,'  arc max ph-bias residual rms      : %.2f m\n',prm.varc.rmsb);
fprintf(f,'  sat. max outlier rate             : %.0f %%\n',prm.vsat.pout*100);
fprintf(f,'  sat. max postfit residual rms     : %.2f m\n',prm.vsat.rmsf);
fprintf(f,'  sat. max ph-bias residual rms     : %.2f m\n',prm.vsat.rmsb);
fprintf(f,'  rcv. max outlier rate             : %.0f %%\n',prm.vrcv.pout*100);
fprintf(f,'  rcv. max postfit residual rms     : %.2f m\n',prm.vrcv.rmsf);
fprintf(f,'  rcv. max ph-bias residual rms     : %.2f m\n',prm.vrcv.rmsb);
fprintf(f,'  sat. max outage of valid obs.     : %.1f sec\n',prm.maxsatout);
fprintf(f,'  rcv. max outage of valid obs.     : %.1f sec\n',prm.maxrcvout);
fprintf(f,'  sat. clock max standard deviation : %.2f ns\n',prm.maxsatcsig);
fprintf(f,'  rcv. clock max standard deviation : %.2f ns\n',prm.maxrcvcsig);
fprintf(f,'  rcv. position max standard dev.   : %.2f m\n',prm.maxrcvpsig);
fprintf(f,'\n');
if any(strcmp(prm.omodel,{'dbldiff','rcvdiff'}))
    fprintf(f,'baselines\n');
    for n=1:size(prm.baseline,1)
        p=readpos(0,0,prm.baseline(n,:),'','approx');
        fprintf(f,'%-7s %-7s : %8.3fkm\n',prm.baseline{n,:},norm(p(1,:,1)-p(1,:,2))/1E3);
    end
    fprintf(f,'\n');
end
fprintf(f,'excluded satellites/receivers\n');
for n=1:size(prm.exclude,1)
    fprintf(f,'%-7s : %s\n',prm.exclude{n,1},pstr(prm.exclude{n,2},prm.exclude{n,3}));
end
fprintf(f,'\n');
satobt={'-','<estimated>','IGS Final','IGS Rapid','IGS Urapid','Est(forwd)','Est(back) ','Est(smooth)','ephs','IGS Pred','COD','EMR','ESA','GFZ','JPL','MIT','COD rapid'};
satclk={'-','<estimated>','IGS Final','IGS Rapid','IGS Urapid','Est(forwd)','Est(back) ','Est(smooth)','IGS Pred','COD','EMR','ESA','GFZ','JPL','MIT','IGS/COD','COD Rapid','IGR/CODR','IGS30s','CODE5s','IGS/CODE5s'};
satsrp={'-','<estimated>'};
rcvpos={'-','<estimated>','IGS00  ','ITRF2000 ','ITRF97 ','GSIPOS','Est(forwd)','Est(back) ','Est(smooth)','IGS00  ','IGb00  ','IGS05','ITRF2005','prmfile'};
rcvclk={'-','<estimated>','IGS Final','IGS Rapid','IGS Urapid','Est(forwd)','Est(back) ','Est(smooth)','JPL    '};
rcvztd={'-','<estimated>','IGS ZPD ','IGS ZPDM','Est(forwd)','Est(back) ','Est(smooth)','MSM Online'};
rcvbcp={'<estimated>','Est(forwd)','Est(back) ','Est(smooth)'};
erpsrc={'-','<estimated>','IGS Final','IERS BulB','IERS BulA','IERS C04',...
        'Est(forwd)','Est(back)','Est(smooth)','CODE','NRCan','ESA','GFZ','JPL','MIT'};
ecosrc={'-','<estimated>'};

fprintf(f,'estimated/fixed parameters\n');
fprintf(f,'  sats  :    orbit        clock      srp param\n');
for n=1:size(prm.satest,1)
    fprintf(f,'  %-6s: %-12s %-12s %-12s\n',prm.satest{n,1},...
            satobt{prm.satest{n,2}+1},satclk{prm.satest{n,4}+1},...
            satsrp{prm.satest{n,3}+1});
end
fprintf(f,'\n');
fprintf(f,'  rcvs  :   position      clock       ztd/grad    phase bias\n');
for n=1:size(prm.rcvest,1)
    fprintf(f,'  %-6s: %-12s %-12s %-12s %-12s\n',prm.rcvest{n,1},...
            rcvpos{prm.rcvest{n,4}+1},rcvclk{prm.rcvest{n,2}+1},...
            rcvztd{prm.rcvest{n,3}+1},rcvbcp{prm.rcvest{n,5}+1});
end
fprintf(f,'\n');
fprintf(f,'  erp   : %-12s  geocent   : %-12s\n',erpsrc{prm.est.erp+1},ecosrc{prm.est.eco+1});
fprintf(f,'\n');
fprintf(f,'measurement noise standard deviations\n');
fprintf(f,'  rcvs  : meas noise(m)\n');
for n=1:size(prm.obssig,1), fprintf(f,'  %-6s: %5.3f\n',prm.obssig{n,1:2}); end
fprintf(f,'\n');
fprintf(f,'apriori standard deviations\n');
fprintf(f,'  sats  : pos/vel(m,m/s) srp1   srp2   srp3   srp4   clock(m,m/s)\n');
for n=1:size(prm.satsig,1)
    fprintf(f,'  %-6s:',prm.satsig{n,1}); fprinte(f,prm.satsig{n,2:end});
end
fprintf(f,'\n');
fprintf(f,'  rcvs  :  clock(m,m/s)  ztd/gradl/gradq(m)   position u/e/n(m)  pbias(m)\n');
for n=1:size(prm.rcvsig,1)
    fprintf(f,'  %-6s:',prm.rcvsig{n,1}); fprinte(f,prm.rcvsig{n,2:end});
end
fprintf(f,'\n');
fprintf(f,'  xp/yp(rad) lod(s) geocent(m) scale(ppb)\n ');
fprinte(f,prm.sig.erp,prm.sig.eco);
fprintf(f,'\n');
fprintf(f,'process noise standard deviations(1/sqrt(sec))\n');
fprintf(f,'  sats  : pos/vel(m,m/s) srp1   srp2   srp3   srp4   clock(m,m/s)\n');
for n=1:size(prm.satprn,1)
    fprintf(f,'  %-6s:',prm.satprn{n,1}); fprinte(f,prm.satprn{n,2:end});
end
fprintf(f,'\n');
fprintf(f,'  rcvs  :  clock(m,m/s)  ztd/gradl/gradq(m)   position u/e/n(m)  pbias(m)\n');
for n=1:size(prm.rcvprn,1)
    fprintf(f,'  %-6s:',prm.rcvprn{n,1}); fprinte(f,prm.rcvprn{n,2:end});
end
fprintf(f,'\n');
fprintf(f,'  xp/yp(rad) lod(s) geocent(m) scale(ppb)\n ');
fprinte(f,prm.prn.erp,prm.prn.eco);
fprintf(f,'\n');

% write observation data statistics --------------------------------------------
function writeedit(f,stats,prm)
if isempty(stats), return, end
tobs=[stats.tobs]; ndat=[stats.ndat]; nobs=[stats.nobs]; nslip=[stats.nslip];
noutl=[stats.noutl]; narc=[stats.narc]; mp1=[stats.mp1]; mp2=[stats.mp2];
loge=[stats.loge]; logs=[stats.logs]; logc=[stats.logc];

fprintf(f,'observation data statistics :\n\n');
fprintf(f,'sat/rcv : total obs ext obs valid obs  slip    arc  mp1(m) mp2(m)      outlier   \n');
for n=1:size(tobs,1)
    nz=sum(ndat(n,:)); nout=sum(noutl(n,:)); mp=[0,0]; r=0;
    if nz>0, mp=sqrt([sum(mp1(n,:)),sum(mp2(n,:))]/nz); r=nout/(nz+nout)*100; end
    fprintf(f,'%-7s : %8d %8d %8d %6d %6d %6.3f %6.3f %6d (%4.1f%%)\n',prm.sats{n},...
            sum(tobs(n,:)),nz,sum(nobs(n,:)),sum(nslip(n,:)),sum(narc(n,:)),mp,nout,r);
end
fprintf(f,'\n');
for n=1:size(tobs,2)
    nz=sum(ndat(:,n)); nout=sum(noutl(:,n)); mp=[0,0]; r=0;
    if nz>0, mp=sqrt([sum(mp1(:,n)),sum(mp2(:,n))]/nz); r=nout/(nz+nout)*100; end
    fprintf(f,'%-7s : %8d %8d %8d %6d %6d %6.3f %6.3f %6d (%4.1f%%)\n',prm.rcvs{n},...
            sum(tobs(:,n)),nz,sum(nobs(:,n)),sum(nslip(:,n)),sum(narc(:,n)),mp,nout,r);
end
fprintf(f,'\n');
nz=sum(sum(ndat)); nout=sum(sum(noutl)); mp=[0,0]; r=0;
if nz>0, mp=sqrt([sum(sum(mp1)),sum(sum(mp2))]/nz); r=nout/(nz+nout)*100; end
fprintf(f,'TOTAL   : %8d %8d %8d %6d %6d %6.3f %6.3f %6d (%4.1f%%)\n',...
        sum(sum(tobs)),nz,sum(sum(nobs)),sum(sum(nslip)),sum(sum(narc)),mp,nout,r);
fprintf(f,'\n');
fprintf(f,'cycle-slip statistics :\n');
if ~isempty(loge)
    fprintf(f,'                                   slips            out    arcs   gfsig  mwsig\n');
    fprintf(f,' sat     rcv   :   obs ion  mw mp1 mp2  if sum fix       org sum   (cm)   (m)\n');
    fprintl(f,loge);
end
fprintf(f,'\ncycle-slip fixed :\n');
if ~isempty(logs)
    fprintf(f,'                                                      slips\n');
    fprintf(f,' sat     rcv            time           n1(cyc)   n2(cyc) wl(cyc) gf(cm) s1   s2\n');
    fprintl(f,logs);
end
fprintf(f,'\nclock-jump fixed :\n');
if ~isempty(logc)
    fprintf(f,'\n');
    fprintf(f,' rcv             time         jump(ms) time-tag offset flg\n');
    fprintl(f,logc);
end
fprintf(f,'\n');

% write receiver info ----------------------------------------------------------
function writeobss(f,stats,prm)
if isempty(stats), return, end
rcvprm=prm_gpsrcvs;
fprintf(f,'receiver names and information :\n\n');
for n=1:length(prm.rcvs)
    i=min(find(strcmp(rcvprm(:,1)',prm.rcvs{n})));
    if ~isempty(i), name=rcvprm{i,3}; info=rcvprm{i,4}; else name=''; info=''; end
    fprintf(f,'%-7s : %-8s %s\n',prm.rcvs{n},name,info);
end
fprintf(f,'\n');
fprintf(f,'receiver/antenna models and antenna deltas h/e/n(m) :\n\n');
for n=1:length(stats)
    fprintf(f,'%-7s : %-22s %-22s %7.4f %7.4f %7.4f\n',stats{n}{:});
end
fprintf(f,'\n');

% write satellite unhealthy/eclipse conditions ---------------------------------
function writeecls(f,stats,prm)
fprintf(f,'satellite unhealthy :\n\n');
for n=1:size(stats,1)
    if (stats(n,5)==1)
        ts=mjdtocal(stats(n,2),stats(n,3)); te=mjdtocal(stats(n,2),stats(n,4));
        fprintf(f,'%-7s : %s-%s\n',prm.sats{stats(n,1)},tstr(ts),tstr(te));
    end
end
fprintf(f,'\n');
fprintf(f,'satellite eclipse :\n\n');
for n=1:size(stats,1)
    if (stats(n,5)==2)
        ts=mjdtocal(stats(n,2),stats(n,3)); te=mjdtocal(stats(n,2),stats(n,4));
        fprintf(f,'%-7s : %s-%s\n',prm.sats{stats(n,1)},tstr(ts),tstr(te));
    end
end
fprintf(f,'\n');

% write estimation statistics --------------------------------------------------
function writeests(f,stats,prm,pass)
if isempty(stats), return, end
ss=zeros(length(prm.sats),7); sr=zeros(length(prm.rcvs),7); st=zeros(1,7);
for n=1:length(stats)
    s=stats{n}{2}; maxo=s(7); s=[s(1:2),s(3:5).^2*s(1),s(6)*s(1)];
    i=find(strcmp(stats{n}{1},prm.sats));
    j=find(strcmp(stats{n}{1},prm.rcvs));
    if ~isempty(i)
        ss(i,1:6)=ss(i,1:6)+s; ss(i,7)=max(ss(i,7),maxo);
        st(1,1:6)=st(1,1:6)+s; st(1,7)=max(st(1,7),maxo);
    elseif ~isempty(j)
        sr(j,1:6)=sr(j,1:6)+s; sr(j,7)=max(sr(j,7),maxo);
    end
end
fprintf(f,'estimation statistics (pass-%d) :\n',pass);
fprintf(f,'                                              residuals(m)            max outage\n');
fprintf(f,'sat/rcv :  total obs outlier prefit-rms postfit-rms ph.b-rms  ph.b-ave   (sec)  \n');
for n=1:length(prm.sats)
    if ss(n,1)>1, ss(n,3:6)=[sqrt(ss(n,3:5)/ss(n,1)),ss(n,6)/ss(n,1)]; end
    fprintf(f,'%-7s : %9d %6d %10.4f %10.4f %10.4f %10.4f %7.0f\n',prm.sats{n},ss(n,:));
end
fprintf(f,'\n');
for n=1:length(prm.rcvs)
    if sr(n,1)>0, sr(n,3:6)=[sqrt(sr(n,3:5)/sr(n,1)),sr(n,6)/sr(n,1)]; end
    fprintf(f,'%-7s : %9d %6d %10.4f %10.4f %10.4f %10.4f %7.0f\n',prm.rcvs{n},sr(n,:));
end
fprintf(f,'\n');
if st(1)>0, st(3:6)=[sqrt(st(3:5)/st(1)),st(6)/st(1)]; end
fprintf(f,'TOTAL   : %9d %6d %10.4f %10.4f %10.4f %10.4f %7.0f\n',st);
fprintf(f,'\n');

% on/off -----------------------------------------------------------------------
function s=onoff(f), if f, s='on'; else s='off'; end
function s=modelclk(f), if f, s='gauss-marcov'; else s='white-noise'; end
function s=modelpos(f), str={'static','kinematic','dynamic','satellite'}; s=str{f+1};
function s=modelatt(f), if f, s='leo satellite'; else s='fixed station'; end
function s=modeltrg(f), str={'none','linear','quadratic'}; s=str{f+1};
function s=modelzpd(f), if ~f, s='random-walk'; else s='gauss-marcov'; end
function s=exepass(p), if ~p.gpsest, s='off'; elseif ~p.backward, s='1-pass';
    elseif ~p.iteration, s='2-pass'; else s='3-pass'; end
function s=obsopt1(f)
    if f==1, s='adjust phase'; elseif f==2, s='adjust code'; else s='off'; end
function s=coropt1(f)
    if f==1, s='sat clock'; elseif f==2, s='satclk+shapiro-delay'; else s='off'; end

% mjd->gpsweek/day/hour --------------------------------------------------------
function wd=gpswd(t)
d=caltomjd(t)-44244; wd(1)=floor(d/7); wd(2)=floor(d-wd(1)*7);

% time string ------------------------------------------------------------------
function s=tstr(ts), s=sprintf('%04d/%02d/%02d %02d:%02d:%02.0f',ts);
function s=pstr(ts,te), s=sprintf('%s-%02d:%02d:%02.0f',tstr(ts),te(4:6));

% prints -----------------------------------------------------------------------
function fprintl(f,log), for n=1:length(log), fprintf(f,'%s\n',log{n}); end

function fprinte(f,varargin)
for v=[varargin{:}]
    s=sprintf('%.1e',v); fprintf(f,' %s%+.0f',s(1:end-4),str2num(s(end-3:end)));
end
fprintf(f,'\n');
