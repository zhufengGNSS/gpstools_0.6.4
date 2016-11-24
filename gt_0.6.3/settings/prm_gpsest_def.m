function prm=prm_gpsest
%-------------------------------------------------------------------------------
% [system] : GpsTools
% [module] : default processing parameters
% [func]   : default processing parameters for gpsedit/gpsestd
% [argin]  : none
% [argout] : prm = processing parameters struct
% [note]   :
% [version]: $Revision: 3 $ $Date: 06/07/14 20:58 $
% [history]: 04/11/20  0.1  new
%            05/03/22  0.2  add parameter prm.tover,prm.sclkmodel,prm.rclkmodel
%                           prm.eclnprn
%-------------------------------------------------------------------------------

% satellites/stations ----------------------------------------------------------
%
%     prm.sats{n,1}   = satellite code ('GPSnn':GPS satellite PRNnn)
%     prm.rcvs{n,1}   = station code
%
%-------------------------------------------------------------------------------
prm.sats={
'GPS01'
'GPS02'
'GPS03'
'GPS04'
'GPS05'
'GPS06'
'GPS07'
'GPS08'
'GPS09'
'GPS10'
'GPS11'
'GPS13'
'GPS14'
'GPS15'
'GPS16'
'GPS17'
'GPS18'
'GPS20'
'GPS21'
'GPS23'
'GPS24'
'GPS25'
'GPS26'
'GPS27'
'GPS28'
'GPS29'
'GPS30'
'GPS31'
};
prm.rcvs={
};

% processing parameters --------------------------------------------------------
%
%     prm.tstart      = observation/estimation start date/time (GPST)
%     prm.tend        = observation/estimation end date/time (GPST)
%     prm.tint        = sampling/estimation interval (sec)
%     prm.tspan       = observation/estimation span (hr)
%     prm.tunit       = processing unit time (hr)
%     prm.tover       = processing overlaping time (hr)
%
%   proc#      tstart                                                 tend
%   -------------+-----------------+-----------------+--------//-------+-------
%          tover |      tunit      | tover           |                 |
%    #1   |<-----+-----------------+----->|          |                 |
%    #2          |          |<-----+-----------------+----->|          |
%    #3          |                            |<-----+--------         |
%    ...         |                                           ....      |
%    #n          |                                              -------+----->|
%                |<- - - - - - - - - - -  tspan - - - - - - - //- - - >|
%
%-------------------------------------------------------------------------------
prm.tstart      = [2004,1,1,0,0,0];
prm.tend        = [2004,1,1,0,0,0];
prm.tint        = 300;
prm.tunit       = 24;
prm.punit       = 0;
prm.tover       = 0;

% output setting ---------------------------------------------------------------
%
%     prm.estout      = output estimation results (1:on)
%     prm.mgeout      = output merged estimation results (1:on)
%     prm.statout     = output statistics/residuals information (1:on)
%     prm.dbout       = output debug information (1:on)
%
%-------------------------------------------------------------------------------
prm.estout      = 1;
prm.mgeout      = 0;
prm.statout     = 0;
prm.dbout       = 0;

% observation data editor settings ---------------------------------------------
%
%     prm.f1          = L1 frequency (GHz)
%     prm.f2          = L2 frequency (GHz)
%     prm.L1code      = L1 code
%                       'C1'   : C1
%                       'P1'   : P1
%                       'P1C1' : P1/C1+(P1-C1)
%     prm.obs.elmin   = min elevation angle (deg) (0:all observation)
%     prm.obs.gapmax  = max arc gap (sec)
%     prm.obs.arcmin  = min arc length (sec)
%     prm.obs.pntmin  = min arc points (ponts)
%     prm.obs.sigmw   = noise of Melbourne-Wubbena (m)
%     prm.obs.sigif   = noise of ion-free code (m)
%     prm.obs.dgmax   = cycle-slip threshold of LG (m)
%     prm.obs.dmp1    = cycle-slip threshold of multilpass L1 (m)
%     prm.obs.dmp2    = cycle-slip threshold of multilpass L2 (m)
%     prm.obs.dwl     = cycle-slip threshold of widelane (cycle)
%     prm.obs.outl    = outlier threshold of MW (sigma)
%     prm.obs.elwe    = threshold elevation angle weighting (0:off,1:on)
%     prm.obs.wind    = cycle-slip detection window length (points)
%     prm.obs.clkrep  = repair clock jump (1:on)
%     prm.obs.reps    = repair cycle-slip (1:on)
%     prm.obs.npnt    = no. of obs. points for fitting
%     prm.obs.nmax    = degree of poly. for fitting
%     prm.obs.separc  = separate arc at eclipse boundary
%     prm.obs.srfilt  = sidereal filter (1:on)
%     prm.obs.srdays  = sidereal filter start/end days relative to current
%     prm.obs.sradj   = sidereal day adjustment parameter (sec)
%     prm.obs.debug   = debug flag (1:on)
%
%-------------------------------------------------------------------------------
prm.f1          = 1.57542;
prm.f2          = 1.22760;
prm.L1code      = 'P1C1';
prm.obs.elmin   = 5;
prm.obs.gapmax  = 300;
prm.obs.arcmin  = 900;
prm.obs.pntmin  = 0;
prm.obs.sigmw   = 1.0;
prm.obs.sigif   = 2.0;
prm.obs.dgmax   = 0.05;
prm.obs.dmp1    = 1;
prm.obs.dmp2    = 1;
prm.obs.dwl     = 0.5;
prm.obs.outl    = 4;
prm.obs.elwe    = 0;
prm.obs.wind    = [30,30];
prm.obs.clkrep  = 1;
prm.obs.reps    = 0;
prm.obs.npnt    = 10;
prm.obs.nmax    = 3;
prm.obs.separc  = 1;
prm.obs.srfilt  = 0;
prm.obs.srdays  = [-1,-1];
prm.obs.sradj   = -9;
prm.obs.debug   = 0;

% parameter estimator settings -------------------------------------------------
%
%     prm.obsedit     = execute observation data editor     (1:on,0:off)
%     prm.gpsest      = execute parameter estimator         (1:on,0:off)
%     prm.loadstate   = use previous est. as initial states (1:on,0:off)
%     prm.filter      = filter
%                          'filterekf' : standard Kalman filter
%                          'filtersrcf': square root conv. filter
%     prm.backward    = backward filter                     (1:on,0:off)
%     prm.iteration   = forward filter iteration            (1:on,0:off)
%     prm.omodel      = estimation strategy
%                          'zerodiff'  : ZD(zero-difference)
%                          'ppp'       : PPP(precise point positioning)
%                          'satdiff'   : SD(single-diff. between satellites.)
%                          'rcvdiff'   : TT(single-diff. between stations)
%                          'dbldiff'   : DD(double difference)
%                          'satzerodiff': SD(satdiff and zerodiff)
%     prm.observ      = Observables
%                          'LC'        : L1/L2 Ionosphere-Free LC
%                          'L1'        : L1 carrier phase
%                          'L2'        : L2 carrier phase
%                          'L1L2'      : L1/L2 carrier phase
%     prm.clkref      = reference clock station             ('':none)
%     prm.tpos        = reference date of station positions
%                       ([] : same as estimation start date)
%     prm.tuinp       = processing unit time (hr) of input states
%
%         back. iter.                  processing passes
%       --------------------------------------------------------------------
%
%         off   off     pass1  +--------------->--------------->
%
%         on    off     pass1  +--------------->--------------->
%                       pass2  <---------------<---------------+
%
%         off   on      pass1  +--------------->--------------->
%                       pass2  +--------------->--------------->
%
%         on    on      pass1  +--------------->--------------->
%                       pass2  <---------------<---------------+
%                       pass3  +--------------->--------------->
%
%-------------------------------------------------------------------------------
prm.obsedit     = 1;
prm.gpsest      = 1;
prm.loadstate   = 0;
prm.filter      = 'filterekf';
prm.backward    = 1;
prm.iteration   = 1;
prm.omodel      = 'zerodiff';
prm.observ      = 'LC';
prm.clkref      = '';
prm.tpos        = [];
prm.batch       = {};
prm.tuinp       = 24;

% baselines --------------------------------------------------------------------
%
%     prm.baseline{n,1}= station1 of baseline for DD/TT
%     prm.baseline{n,2}= station2 of baseline for DD/TT ('ALL': all stations)
%
%-------------------------------------------------------------------------------
prm.baseline={
};

% excluded satellites/stations -------------------------------------------------
%
%     prm.exclude{n,1}= excluded satellite/station
%     prm.exclude{n,2}= excluded time start [year,month,day,hour,min,sec]
%     prm.exclude{n,3}= excluded time end   [year,month,day,hour,min,sec]
%
%-------------------------------------------------------------------------------
prm.exclude={
};

% input/output data settings ---------------------------------------------------
%
%     prm.dirs.obs    = observation data directory
%     prm.dirs.obc    = pre-processed observation data directory
%     prm.dirs.nav    = navigation message directory
%     prm.dirs.eph    = satellite ephemeris directory
%     prm.dirs.clk    = satellite/station clock directory
%     prm.dirs.pos    = staion position directory
%     prm.dirs.erp    = earth rotation parameters directory
%     prm.dirs.trop   = tropospheric parameters directory
%     prm.dirs.ion    = ionospheric parameters directory
%     prm.dirs.inp    = input satellite estimation results directory
%     prm.dirs.inpr   = input receiver estimation results directory
%     prm.dirs.est    = output estimation results directory
%
%         following keywords replaced in directory paths
%
%               %P -> install path
%               %S -> satellite/station name (UPPER case)
%               %s -> satellite/station name (lower case)
%               %r -> station name (only valid for prm.dirs.est)
%               %G -> GSI station code (tail 4 characters)
%               %Y -> year           (1960-2099)
%               %y -> year           (00-99)
%               %m -> month          (01-12)
%               %d -> day            (01-31)
%               %n -> day of year    (001-366)
%               %W -> gpsweek no     (0000-9999)
%               %D -> day of gpsweek (0-6)
%               %M -> iers bulletin month (1- )
%
%     prm.obs.src     = input observation data
%                           'rinex'   : RINEX OBS daily (24hr)
%                           'rinex3'  : RINEX OBS 3 hours (3hr)
%                           'rinex1'  : RINEX OBS hourly (1hr)
%                           'rinexh'  : RINEX OBS high-rate (15min)
%                           'obsg'    : simulated obs. data
%     prm.src.nav     = input navigation message
%                           'brdc'    : IGS combined RINEX NAV (24hr)
%                           'rinex'   : RINEX NAV daily (24hr)
%                           'rinex3'  : RINEX NAV 3 hours (3hr)
%                           'rinex1'  : RINEX NAV hourly (1hr)
%                           'rinexh'  : RINEX NAV high-rate (15min)
%-------------------------------------------------------------------------------
prm.dirs.obs    = '';
prm.dirs.obc    = '';
prm.dirs.nav    = '';
prm.dirs.eph    = '';
prm.dirs.clk    = '';
prm.dirs.pos    = '';
prm.dirs.erp    = '';
prm.dirs.trop   = '';
prm.dirs.ion    = '';
prm.dirs.inp    = '';
prm.dirs.inpr   = '';
prm.dirs.est    = '';
prm.obs.src     = 'rinex';
prm.src.nav     = 'rinex';

% estimation/measurement model settings ----------------------------------------
%
%     prm.elmin       = min elevation angle (deg)
%     prm.elmax       = max elevation angle (deg)
%     prm.trop        = tropospheric delay model
%                          'trop_saast': Saastamoinen
%                          'trop_mhop' : modefied Hopfield
%     prm.mapf        = tropospheric mapping function
%                          'mapf_cosz' : cos z
%                          'mapf_nmf'  : Niell mapping function
%                          'mapf_gmf'  : Global mapping function
%     prm.zpdmodel    = tropospheric zenith path delay model
%                          0           : random-walk
%                          1           : 1st order gauss-marcov
%     prm.trgmodel    = tropospheric gradient model
%                          0           : none
%                          1           : linear gradient
%                          2           : quadratic gradient
%     prm.sclkmodel   = satellite clock model
%                          0           : white-noise
%                          1           : 1st order gauss-marcov
%     prm.rclkmodel   = receiver clock model
%                          0           : white-noise
%                          1           : 1st order gauss-marcov
%     prm.rposmodel   = receiver position model
%                          0           : static model
%                          1           : kinematic model
%                          2           : dynamic model
%     prm.rattmodel   = receiver attitude model
%                          0           : fixed station
%                          1           : north:alongtrk,west:crosstrk,up:radial
%     prm.sigweight   = meas. noise weighting by elevation angle (1:on,0:off)
%     prm.eclnsig     = meas. noise factor in eclipse & post-eclipse maneuver
%                         prm.eclnsig(1) : Block II/IIA
%                         prm.eclnsig(2) : Block IIR
%     prm.eclnprn     = process noise factor in post-eclipse maneuver
%                         prm.eclnprn(1) : Block II/IIA
%                         prm.eclnprn(2) : Block IIR
%     prm.ecltime     = max post-eclipse maneuver time (sec)
%     prm.ogapmax     = reset time of observation outage (sec)
%     prm.chkcov      = cov. positive definiteness test (1:on,0:off)
%     prm.maxiter     = max count of filter iteration
%     prm.minobs      = min counts of observation data for estimation (0:any)
%     prm.sdfact      = std. dev. factor of initial states
%                       (0:reset to a priori std. dev)
%     prm.sdfact2     = std. dev. factor of pre-determined states
%                       (0:do not use std. dev/accuracy code)
%     prm.clkcorr     = correct clock by phase bias residual (1:on,0:off)
%     prm.clkqc       = qc clock compared with reference clock
%                          ''        : off
%                          'igs'     : IGS Final
%                          'igr'     : IGR Rapid
%                          'code'    : CODE
%     prm.clkalign    = align clock to reference clock
%                          ''        : off
%                          'igs'     : IGS Final
%                          'igr'     : IGR Rapid
%     prm.corrf       = mesurement corrections (1:on,0:off)
%                          prm.corrf(1) : satellite antenna offsets
%                          prm.corrf(2) : receiver antenna offsets
%                          prm.corrf(3) : relativistic effects
%                          prm.corrf(4) : phase windup
%                          prm.corrf(5) : phase multipath
%     prm.sitedisp    = station displacement corrections (1:on,0:off)
%                          prm.sitedisp(1) : solid earth tides
%                          prm.sitedisp(2) : ocean loading
%                          prm.sitedisp(3) : polar tides
%                          prm.sitedisp(4) : eliminate permanent deformation
%     prm.erpvar      = erp variation correction (1:on,0:off)
%     prm.utc_tai     = utc - tai (sec)
%     prm.nutmodel    = nutation model
%                          ''        : iau1976/iau1980
%                          'iers1996': iers1996
%     prm.metsrc      = meteorological parameters
%                          ''          : use default value
%                          'mso'       : read msm gpv data
%     prm.metprm      = default meteorological parameters
%                          prm.metprm(1) : pressure (hPa)
%                          prm.metprm(2) : temperature (C)
%                          prm.metprm(3) : retive humidity (%)
%     prm.initerp     = initial earth rotation parameters
%                          'igs'  : IGS Final
%                          'bulb' : IERS bulletin B
%                          'c04'  : IERS C04
%                          'erpf' : estimated (forward)
%                          'erpb' : estimated (backward)
%     prm.clkintp     = clock interpolation (1:on,0:off)
%
%-------------------------------------------------------------------------------
prm.elmin       = 10;
prm.elmax       = 90;
prm.trop        = 'trop_saast';
prm.mapf        = 'mapf_nmf';
prm.zpdmodel    = 1;
prm.trgmodel    = 0;
prm.sclkmodel   = 0;
prm.rclkmodel   = 0;
prm.rposmodel   = 0;
prm.rattmodel   = 0;
prm.sigweight   = 0.6;
prm.eclnsig     = [0,1];
prm.eclnprn     = [10,10];
prm.ecltime     = 1800;
prm.ogapmax     = 7200;
prm.chkcov      = 0;
prm.maxiter     = 1;
prm.minobs      = 0;
prm.sdfact      = 1;
prm.sdfact2     = 0;
prm.clkcorr     = 1;
prm.clkqc       = '';
prm.clkalign    = '';
prm.corrf       = [1,1,1,1,0];
prm.sitedisp    = [1,1,1,0];
prm.erpvar      = 1;
prm.utc_tai     = -32;
prm.nutmodel    = 'iers1996';
prm.metsrc      = '';
prm.metprm      = [1013.25,15,50];
prm.initerp     = '';
prm.clkintp     = 0;

% satellite orbit model setting ------------------------------------------------
%
%     prm.obt.g_nmax  = max degree of geopotential
%     prm.obt.p_plgrv = solar/planetary potential
%                           ''          : none
%                           'grv_sunmoon': sun and moon
%     prm.obt.p_solarpr=solar radiation pressure model
%                           ''          : none
%                           'srp_simple': simple model
%                           'srp_rock4' : ROCK4/42 model
%                           'srp_gspm'  : GSPM model
%                           'srp_gspmm' : GSPM+Along/Cross model
%                           'srp_code'  : CODE model(3 parameters)
%                           'srp_code2' : CODE model(4 parameters)
%                           'srp_code3' : CODE model(6 parameters)
%     prm.obt.p_eclipse=eclipse model
%                           ''          : cylindric model
%                           'penumbra'  : penumbra/umbra model
%     prm.obt.p_atmos = atmospheric drag model
%                           ''          : none
%     prm.obt.p_tidal = earth tides
%                           ''          : none
%                           'solid'     : solid earth tides
%                           'otide'     : solid eath and ocean tides
%     prm.obt.p_relativ=relativistic effects
%                           ''          : none
%                           'iers'      : iers conventions
%     prm.obt.p_deltav= thrust forces/delta-v
%                           ''          : none
%                           'dv_on'     : thruster on
%     prm.obt.modelnut= precession/nutation model
%                           ''          : iau1976/iau1980
%     prm.obt.tstep   = integration step time (sec)
%
%-------------------------------------------------------------------------------
prm.obt.g_nmax   = 8;
prm.obt.p_plgrv  = 'grv_sunmoon';
prm.obt.p_solarpr= 'srp_code';
prm.obt.p_eclipse= '';
prm.obt.p_atmos  = '';
prm.obt.p_tidal  = 'solid';
prm.obt.p_relativ= 'iers';
prm.obt.p_deltav = '';
prm.obt.tstep    = 30;
prm.obt.modelnut = '';

% quality control settings -----------------------------------------------------
%
%     prm.reests      = re-estimate without excluded satellites   (1:on,0:off)
%     prm.reestr      = re-estimate without excluded receivers    (1:on,0:off)
%     prm.outlp       = outlier threshold of prefit res. (sigma)  (0:no check)
%     prm.outlf       = outlier threshold of postfit res. (sigma) (0:no check)
%     prm.varc.pout   = arc max outlier rate                      (0:no check)
%     prm.varc.rmsf   = arc max postfit residual rms (m)          (0:no check)
%     prm.varc.rmsb   = arc max phase bias residual rms (m)       (0:no check)
%     prm.vsat.pout   = sat. max outlier rate                     (0:no check)
%     prm.vsat.rmsf   = sat. max postfit residual rms (m)         (0:no check)
%     prm.vsat.rmsb   = sat. max phase bias residual rms (m)      (0:no check)
%     prm.vrcv.pout   = rcv. max outlier rate                     (0:no check)
%     prm.vrcv.rmsf   = rcv. max postfit residual rms (m)         (0:no check)
%     prm.vrcv.rmsb   = rcv. max phase bias residual rms (m)      (0:no check)
%     prm.maxsatout   = sat. max outage of valid obs (sec)        (0:no check)
%     prm.maxrcvout   = rcv. max outage of valid obs (sec)        (0:no check)
%     prm.maxsatcsig  = max satellite clock sigma (ns)            (0:no check)
%     prm.maxrcvcsig  = max receiver clock sigma (ns)             (0:no check)
%     prm.maxrcvpsig  = max receiver position sigma (m)           (0:no check)
%
%-------------------------------------------------------------------------------
prm.reests      = 1;
prm.reestr      = 1;
prm.outlp       = 4;
prm.outlf       = 4;
prm.varc.pout   = 0;
prm.varc.rmsf   = 0;
prm.varc.rmsb   = 0;
prm.vsat.pout   = 0;
prm.vsat.rmsf   = 0;
prm.vsat.rmsb   = 0;
prm.vrcv.pout   = 0;
prm.vrcv.rmsf   = 0;
prm.vrcv.rmsb   = 0;
prm.maxsatout   = 0;
prm.maxrcvout   = 0;
prm.maxsatcsig  = 0;
prm.maxrcvcsig  = 0;
prm.maxrcvpsig  = 0;

% external parameter files -----------------------------------------------------
%
%     prm.satsrpf     = satellite srp parameters file path (m-file)
%     prm.rcvposf     = station position file path (m-file)
%     prm.pcv         = receiver antenna pcv parameters file path
%     prm.satpcv      = satellite antenna pcv parameters file path
%     prm.oload       = ocean loading parameters file path
%     prm.mpc         = code multipath profile (mat-file)
%     prm.mpp         = phase multipath profile (mat-file)
%     prm.exsatrcv    = excluded satellites/receivers on estimation (m-file)
%
%-------------------------------------------------------------------------------
prm.satsrpf     = '';
prm.rcvposf     = '';
prm.pcv         = '%P\data\igs_01.pcv';
prm.satpcv      = '';
prm.oload       = '%P\data\oload_igs.blq';
prm.mpc         = '';
prm.mpp         = '';
prm.exsatrcv    = '';

% estimated/fixed parameters ---------------------------------------------------
%
% satellite parameters
%     prm.satest{n,1} = satellite ('ALL':all satellites)
%     prm.satest{n,2} = satellite orbit
%                          0 : fixed to broadcast ephemeris
%                          1 : <estimated>
%                          2 : fixed to IGS final
%                          3 : fixed to IGS rapid
%                          4 : fixed to IGS ultra-rapid
%                          5 : fixed to estimated (forward) by gpsestd
%                          6 : fixed to estimated (backward) by gpsestd
%                          7 : fixed to estimated (smoothed) by gpsestd
%                          8 : fixed to ephemeris generated by geneph
%                          8 : fixed to IGS ultra-rapid (predicted)
%                         10 : fixed to CODE
%                         11 : fixed to NRCan
%                         12 : fixed to ESA
%                         13 : fixed to GFZ
%                         14 : fixed to JPL
%                         15 : fixed to MIT
%                         16 : fixed to CODE rapid
%     prm.satest{n,3} = satellite srp parameters
%                          0 : fixed to satellite parameters
%                          1 : <estimated>
%     prm.satest{n,4} = satellite clock
%                          0 : fixed to broadcast ephemeris
%                          1 : <estimated>
%                          2 : fixed to IGS final
%                          3 : fixed to IGS rapid
%                          4 : fixed to IGS ultra-rapid
%                          5 : fixed to estimated (forward) by gpsestd
%                          6 : fixed to estimated (backward) by gpsestd
%                          7 : fixed to estimated (smoothed) by gpsestd
%                          8 : fixed to IGS ultra-rapid (predicted)
%                          9 : fixed to CODE
%                         10 : fixed to NRCan
%                         11 : fixed to ESA
%                         12 : fixed to GFZ
%                         13 : fixed to JPL
%                         14 : fixed to MIT
%                         15 : fixed to IGS/COD (IGS final+CODE interpolation)
%                         16 : fixed to CODE rapid
%                         17 : fixed to IGR/CODR (IGS rapid+CODE rapid interpolation)
% station parameters
%     prm.rcvest{n,1} = station ('ALL':all stations)
%     prm.rcvest{n,2} = receiver clock
%                          0 : fixed to point positioning results
%                          1 : <estimated>
%                          2 : fixed to IGS final
%                          3 : fixed to IGS rapid
%                          4 : fixed to IGS ultra-rapid
%                          5 : fixed to estimated (forward) by gpsestd
%                          6 : fixed to estimated (backward) by gpsestd
%                          7 : fixed to estimated (smoothed) by gpsestd
%                          8 : fixed to IGS ultra-rapid (predicted)
%                          9 : fixed to CODE
%                         10 : fixed to NRCan
%                         11 : fixed to ESA
%                         12 : fixed to GFZ
%                         13 : fixed to JPL
%                         14 : fixed to MIT
%     prm.rcvest{n,3} = tropospheric zenith total delay/gradients
%                          0 : fixed to tropospheric model
%                          1 : <estimated>
%                          2 : fixed to igs zpd
%                          3 : fixed to igs zpd monthly
%                          4 : fixed to estimated (forward) by gpsestd
%                          5 : fixed to estimated (backward) by gpsestd
%                          6 : fixed to estimated (smoothed) by gpsestd
%                          7 : fixed to gpv derived zpd (msm online)
%     prm.rcvest{n,4} = receiver position
%                          0 : fixed to point positioning results
%                          1 : <estimated>
%                          2 : fixed to IGS final
%                          3 : fixed to ITRF2000
%                          4 : fixed to ITRF97
%                          5 : fixed to GSI estimated position
%                          6 : fixed to estimated (forward) by gpsestd
%                          7 : fixed to estimated (backward) by gpsestd
%                          8 : fixed to estimated (smoothed) by gpsestd
%                          9 : fixed to IGS00
%                         10 : fixed to IGb00
%                         11 : fixed to position read from parameter file
%     prm.rcvest{n,5} = phase bias
%                          0 : <estimated>
%                          1 : fixed to estimated (forward)
%                          2 : fixed to estimated (backward)
%                          3 : fixed to estimated (smoothed)
%
% global parameters
%     prm.est.erp     = earth rotation parameters
%                          0 : fixed to initial ERP
%                          1 : <estimated>
%                          2 : fixed to IGS Final
%                          3 : fixed to IGS rapid
%                          4 : fixed to IERS bulletin B
%                          5 : fixed to IERS bulletin A
%                          6 : fixed to IERS C04
%                          7 : fixed to estimated (forward) by gpsestd
%                          8 : fixed to estimated (backward) by gpsestd
%                          9 : fixed to estimated (smoothed) by gpsestd
%                         10 : fixed to CODE
%                         11 : fixed to NRCan
%                         12 : fixed to ESA
%                         13 : fixed to GFZ
%                         14 : fixed to JPL
%                         15 : fixed to MIT
%     prm.est.eco     = earth mass center offset
%                          0 : fixed to zero
%                          1 : <estimated>
%
%-------------------------------------------------------------------------------
prm.satest={
'ALL',      2, 0, 2
};
prm.rcvest={
'ALL',      1, 1, 1, 1
};
prm.est.erp = 2;
prm.est.eco = 0;
prm.bcpest  = 2;

% a-priori std. deviations -----------------------------------------------------
%
% satellite states
%     prm.satsig{n,1}   = satellite code
%     prm.satsig{n,2}   = a-priori std. dev. of sat. position(m)
%     prm.satsig{n,3}   = a-priori std. dev. of sat. velocity(m/s)
%     prm.satsig{n,4:7} = a-priori std. dev. of srp paramters
%     prm.satsig{n,8:9} = a-priori std. dev. of sat. clock(m,m/s) [bias,drift]
%
% station states
%     prm.rcvsig{n,1}   = station code
%     prm.rcvsig{n,2:3} = a-priori std. dev. of rcv. clock(m,m/s) [bias,drift]
%     prm.rcvsig{n,4:5} = a-priori std. dev. of tropos zpd (m,m/h) [zpd,drift]
%     prm.rcvsig{n,6}   = a-priori std. dev. of tropos gradient (m)
%     prm.rcvsig{n,7:9} = a-priori std. dev. of sta. position(m) [up,east,north]
%     prm.rcvsig{n,10}  = a-priori std. dev. of phase bias(m)
%
% earth rotation parameters
%    prm.sig.erp(1,1)   = a-priori std. dev. of pole offset xp, yp (rad)
%    prm.sig.erp(1,2)   = a-priori std. dev. of ut1-utc(sec)
%
% earth mass center offset
%    prm.sig.eco(1,1)   = a-priori std. dev. of geocenter offset (m)
%    prm.sig.eco(1,2)   = a-priori std. dev. of geocenter scale (ppb)
%
%-------------------------------------------------------------------------------
prm.satsig={
'ALL',   1E+1, 5E-4, 5E-1, 5E-1, 5E-1, 5E-1, 5E+0, 5E-4
};
prm.rcvsig={
'ALL',   1E+1, 1E-3, 3E-1, 1E-2, 1E-2, 0, 0, 0, 1E+0
};
prm.sig.erp=[1E-8, 2E-4];
prm.sig.eco=[1E-2, 5];

% process noises (1/sqrt(sec)) -------------------------------------------------
%
% satellite states
%
%     prm.satprn{n,1}   = satellite code ('ALL':all satellites)
%     prm.satprn{n,2}   = process noise of sat. position(m) [x,y,z]
%     prm.satprn{n,3}   = process noise of sat. velocity(m/s) [x,y,z]
%     prm.satprn{n,4:7} = process noise of srp paramters
%     prm.satprn{n,8:9} = process noise of sat. clock(m,m/h) [bias,drift]
%
% station states
%
%     prm.rcvprn{n,1}   = station code ('ALL':all stations)
%     prm.rcvprn{n,2:3} = process noise of rcv. clock(m,m/h) [bias,drift]
%     prm.rcvprn{n,4:5} = process noise of tropos zpd (m,m/h) [zpd,drift]
%     prm.rcvprn{n,6}   = process noise of tropos gradient (m)
%     prm.rcvprn{n,7:9} = process noise of sta. position(m) [up,east,north]
%
%     prm.prn.ecl(1,1)  = process noise of eclipsing sat. position(m)
%     prm.prn.ecl(1,2)  = process noise of eclipsing sat. velocity(m/s)
%
% earth rotation parameters
%
%    prm.prn.erp(1,1)   = process noise of pole offset xp, yp (rad)
%    prm.prn.erp(1,2)   = process noise of ut1-utc (sec)
%
% earth mass center offset
%
%    prm.prn.eco(1,1)   = process noise of geocenter offset (m)
%    prm.sig.eco(1,2)   = process noise of geocenter scale (ppb)
%
%-------------------------------------------------------------------------------
prm.satprn={
'ALL',   0E+0, 1E-9, 1E-6, 1E-6, 1E-6, 1E-6, 1E-2, 1E-7
};
prm.rcvprn={
'ALL',   1E+1, 1E-2, 1E-4, 1E-5, 1E-5, 0E+0, 0E+0, 0E+0, 1E-4
};
prm.prn.ecl=[0E+0, 5E-9];
prm.prn.erp=[1E-11, 2E-7];
prm.prn.eco=[1E-6, 1E-3];

% measurement noise std. deviations ---------------------------------------------
%
%     prm.obssig{n,1}   = station code ('ALL':all stations)
%     prm.obssig{n,2}   = std. dev. of phase measurement noise(m)(zenith)
%
%-------------------------------------------------------------------------------
prm.obssig={
'ALL',   0.010
};
