%  GpsTools 
%
%  GPS/GNSS precise analysis software package. The package estimates receiver
%  kinematic/static positions, satellite orbit, satellite/receiver clock,
%  tropospheric parameters, using dual-frequency code and carrier phase
%  observables. Estimator employes extended kalman filter/smoother to support
%  analysis strategy of ZD/PPP and DD. Includes gui for seting parameters,
%  executing analysis, downloader for data/products, and ploting results. 
%
%            Copyright(c) 2004-2008 by T.Takasu, all rights reserved
%
% commands/libraries
%
%   gpstools          (m)     : main menu
%   gpsfiles          (m)     : file manipulation tool
%   ggpsest           (m)     : parameter estimator/obs data editor gui
%   gpsestd           (m)     : parameter estimator
%   gpsdown           (m)     : data/products downloader
%   plotobs           (m)     : plot observation data
%   plotpos           (m)     : plot receiver positions
%   ploteph           (m)     : plot satellite orbits
%   plotclk           (m)     : plot satellite/receiver clocks
%   plottrop          (m)     : plot tropospheric parameters
%   plotbcp           (m)     : plot phase biases
%   plotgpv           (m)     : plot gpv data
%   plotion           (m)     : plot ionospheric parameters
%   plotleo           (m)     : plot leo satellite orbits
%   plotstats         (m)     : plot residuals
%   plotmp            (m)     : plot multipath profile
%   plotpmap          (m)     : plot pwv parameters
%   gpstrack          (m)     : plot satellite ground tracks and receiver positions
%   genblq            (m)     : generate ocean loading parameters
%   geneph            (m)     : interpolate and generate satellite ephemerides
%   pwvgps            (m)     : generate gps pwv
%   pwvzonde          (m)     : generate zonde pwv
%   ephtosp3          (m)     : covert estimation results to sp3
%
%   readobs           (m)     : read observation data
%   readnav           (m)     : read navigation messages
%   readeph           (m)     : read satellite ephemerides
%   readclk           (m)     : read satellite/receiver clocks
%   readpos           (m)     : read receiver positions
%   readtrop          (m)     : read tropospheric parameters
%   readerp           (m)     : read earth rotation parameters
%   readstats         (m)     : read residuals
%   readest           (m)     : read estimation results
%   readpcv           (m)     : read satellite/receiver antenna parameters
%   readmet           (m)     : read meterological parameters
%   readgpv           (m)     : read gpv data
%   readdcb           (m)     : read dcb data
%   readrcv           (m)     : read station information
%   readoload         (m)     : read ocean loading parameters
%
%   readsp3           (mex)   : read sp3 ephemeris/clock file
%   readrinexobs      (mex)   : read rinex obs file
%   readrinexnav      (mex)   : read rinex nav file
%   readrinexclk      (mex)   : read rinex clk file
%   readsinexpos      (m)     : read sinex pos file
%   readionex         (m)     : read ionex file
%   readgsipos        (m)     : read GSI position esitmation file
%   readantex         (m)     : read antex antenna parameter file
%   readgrib          (mex)   : read grib gpv data file
%   readgrace1b       (m)     : read GRACE Level 1B data file
%   readchorb         (m)     : read CHORB CHAMP orbit file
%   readvmfgrid       (mex)   : read VMF grid file
%   readvmf           (m)     : read VMF coefficients
% 
%   satrange          (mex)   : range of satellite-receiver
%   satazel           (mex)   : satellite direction from station
%   isrange           (mex)   : range of inter-satellite
%   navtostate        (mex)   : sat position/velocity from navigation messages
%   pointpos          (mex)   : determin position and clock by point positioning
%   sunmoonpos        (m)     : get solar/lunar positions
%
%   rcvapc            (mex)   : receiver antenna offset correction
%   rcvmpc            (mex)   : multipath bias corrections
%   satapc            (mex)   : satellite antenna offset correction
%   relcorr           (mex)   : relativistic effect correction
%   sitedisp          (mex)   : station displacement correction
%   phwindup          (mex)   : phase-windup correction
%   erpvar            (m)     : erp variation correction
%   esthelmert        (m)     : estimate helmert transformation parameters
%   interpeph         (m)     : interpolate satellite ephemeris
%   bpfilt            (m)     : bandpass filter
%   srfilt            (m)     : sidereal filter
%
%   state_2body       (mex)   : satellite orbit model by 2-body problem
%   state_kepler      (mex)   : satellite orbit model by Keplar
%   state_geoj2       (mex)   : satellite orbit model by J2 precession
%   state_precorbit   (m)     : satellite orbit model by presice purturbation
%   state_satclock    (mex)   : satellite clock model
%
%   trop_saast        (mex)   : saastamoinen tropospheric model
%   trop_mhop         (m)     : modefied hopfield tropospheric model
%   trop_gpt          (mex)   : global pressure and temperature model
%   mapf_cosz         (m)     : cosz mapping function
%   mapf_nmf          (mex)   : niell mapping function (NMF)
%   mapf_gmf          (mex)   : global mapping function (GMF)
%   mapf_vmf1         (mex)   : vienna mapping function 1 (VMF1)
%
%   ion_klob          (mex)   : klobuchar ionosphere model
%   ion_tec           (m)     : tec map ionosphere model
%
%   caltomjd          (mex)   : calender date/time to mjd
%   mjdtocal          (mex)   : mjd to calender date/time
%   ecsftoecef        (mex)   : transform eci position to ecef
%   geodtoecef        (mex)   : transform geodetic position to ecef
%   eceftogeod        (mex)   : transform ecef position to geodetic
%   ecsftosatf        (mex)   : transform eci to satellite orbit coord.
%   geomtogeoc        (mex)   : geomagnetic to geocentric postion
%   eletostate        (mex)   : satellite orbit element to states
%   statetoele        (mex)   : satellite state to orbit element
%   sphfunc           (mex)   : spheric harmonic function
%   shadowfunc        (mex)   : shadow function of satellite
%   helmert           (m)     : helmert transformation
%   interplag         (m)     : lagrange interpolation
% 
%   gut               (m)     : gui common libraries for gpstools
%   gmt               (m)     : map common libraries for gpstools
%   ggt               (m)     : graph common libraries for gpstools
%
% [history]: 06/07/23  0.6.3  new
%            08/11/30  0.6.4  add trop_gpt readvmfgrid readvmf mapf_vmf1 (gt_0.6.4)
%                             delete lambda
%
