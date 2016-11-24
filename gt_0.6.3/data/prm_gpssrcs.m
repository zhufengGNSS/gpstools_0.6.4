function addrs=prm_gpssrcs
%------------------------------------------------------------------------------
% [system] : GpsTools
% [module] : GPS data sources definition
% [func]   : GPS data sources definition
% [argin]  : none
% [argout] : addrs : data sources
%                addrs{n,1} : data category(1:OBS,2:PRODUCTS,3:EOP,4:OTHERS)
%                addrs{n,2} : data label
%                addrs{n,3} : data source address
% [note]   : keywords in address replaced as follows.
%                %S -> station/centre code (UPPER case)
%                %s -> station/centre code (lower case)
%                %G -> GSI station code
%                %Y -> year        (1990-2010)
%                %y -> year        (00-99)
%                %m -> month       (01-12)
%                %d -> day         (01-31)
%                %n -> day of year (001-366)
%                %h -> hour        (00-23)
%                %ha-> 3H hour     (00,03,06,...,21)
%                %hb-> 6H hour     (00,06,12,18)
%                %hc-> 12H hour    (00,12)
%                %H -> hour code   (a,b,c,...,x)
%                %t -> minute      (00,15,30,45)
%                %g -> 3H hour     (02,05,08,...,23)
%                %f -> 3H sequence no (1=00,2=03,...,8=21)
%                %W -> gpsweek no. (0000-9999)
%                %N -> sequence no. (0000-9999)
%                %D -> day of gpsweek (0-6)
%                %M -> iers bulletin month (1- )
%                %RO-> RINEX OBS (%s%n0.%yo)
%                %RD-> RINEX COBS (%s%n0.%yd)
%                %RN-> RINEX NAV (GPS) (%s%n0.%yn)
%                %RG-> RINEX NAV (GLONASS) (%s%n0.%yg)
%                %RM-> RINEX MET (%s%n0.%ym)
%                %RS-> RINEX SUM (%s%n0.%ys)
%                %{s1|s2|...|sn}->strings s1,s2,...,sn
% [version]: $Revision: 1 $ $Date: 06/07/08 9:21 $
% [history]: 04/11/18  0.1  new
%            05/03/30  0.2  change format
%            05/04/21  0.3  change gsi directory
%            05/08/11  0.4  add GRACE Level1B
%            06/07/01  0.5  add CODE products
%------------------------------------------------------------------------------

% server address --------------------------------------------------------------
cddis='ftp://cddis.gsfc.nasa.gov';       % CDDIS
sopac='ftp://garner.ucsd.edu/pub';       % SOPAC
ign  ='ftp://igs.ensg.ign.fr/pub';       % IGN
gsi  ='ftp://terras.gsi.go.jp/data';     % GSI
usno ='ftp://tycho.usno.navy.mil/pub';   % USNO
gfd  ='ftp://dennou-k.gfd-dennou.org';   % GFD_Dennou_Club
jma  ='ftp://ddb.kishou.go.jp/pub';      % JMA DDB
podaac='ftp://podaac.jpl.nasa.gov/pub';  % JPL PODACC
code ='ftp://ftp.unibe.ch/aiub/CODE';    % CODE

addrs={...

% observation data ------------------------------------------------------------

1,'IGS OBS DAILY (CDDIS)',  [cddis,'/gps/data/daily/%Y/%n/%yd/%RD.Z']
1,'IGS NAV COMB  (CDDIS)',  [cddis,'/gps/data/daily/%Y/%n/%yn/brdc%n0.%yn.Z']
1,'IGS NAV DAILY (CDDIS)',  [cddis,'/gps/data/daily/%Y/%n/%yn/%RN.Z']
1,'IGS MET DAILY (CDDIS)',  [cddis,'/gps/data/daily/%Y/%n/%ym/%RM.Z']

1,'IGS OBS HOURLY (CDDIS)', [cddis,'/gps/data/hourly/%Y/%n/%h/%s%n%H.%yd.Z']
1,'IGS NAV HOURLY (CDDIS)', [cddis,'/gps/data/hourly/%Y/%n/%h/%s%n%H.%yn.Z']
1,'IGS MET HOURLY (CDDIS)', [cddis,'/gps/data/hourly/%Y/%n/%h/%s%n%H.%ym.Z']

1,'IGS OBS H-RATE (CDDIS)', [cddis,'/gps/data/highrate/%Y/%n/%yd/%h/%s%n%H%t.%yd.Z']
1,'IGS NAV H-RATE (CDDIS)', [cddis,'/gps/data/highrate/%Y/%n/%yd/%h/%s%n%H%t.%yn.Z']
1,'IGS MET H-RATE (CDDIS)', [cddis,'/gps/data/highrate/%Y/%n/%h/%yd/%s%n%H%t.%ym.Z']

1,'IGS OBS DAILY (SOPAC)',  [sopac,'/rinex/%Y/%n/%RD.Z']
1,'IGS MET DAILY (SOPAC)',  [sopac,'/met/%Y/%n/%RM.Z']

1,'GSI OBS 24H 電子基準点', [gsi,'/GPS_products/%Y/%n/%G%n0.%yo.gz']
1,'GSI NAV 24H 電子基準点', [gsi,'/GPS_products/%Y/%n/%G%n0.%yn.gz']
1,'GSI OBS  3H 電子基準点', [gsi,'/GPS_products/%Y/%n/%G%n%f.%yo.gz']
1,'GSI NAV  3H 電子基準点', [gsi,'/GPS_products/%Y/%n/%G%n%f.%yn.gz']

1,'CHAMP OBS (CDDIS)',      [cddis,'/gps/data/satellite/champ/%Y/%n/cham%n0.%yd.Z']
1,'CHAMP SUM (CDDIS)',      [cddis,'/gps/data/satellite/champ/%Y/%n/cham%n0.%ys.Z']
1,'JASON-1 OBS (CDDIS)',    [cddis,'/gps/data/satellite/champ/%Y/%n/jas1%n0.%yd.Z']
1,'JASON-1 SUM (CDDIS)',    [cddis,'/gps/data/satellite/champ/%Y/%n/jas1%n0.%ys.Z']

1,'USNO NANUs',             [usno,'/gps/gpsnanu.%Y']
1,'USNO NANUs sat index',   [usno,'/gps/gpsout%y.txt']

% products ---------------------------------------------------------------------

2,'IGS EPH Final (CDDIS)',  [cddis,'/gps/products/%W/igs%W%D.sp3.Z']
2,'IGS CLK Final (CDDIS)',  [cddis,'/gps/products/%W/igs%W%D.clk.Z']
2,'IGS ERP Final (CDDIS)',  [cddis,'/gps/products/%W/igs%W7.erp.Z']
2,'IGS SUM Final (CDDIS)',  [cddis,'/gps/products/%W/igs%W7.sum.Z']
2,'IGS EPH Final(SP3C)(CDDIS)', [cddis,'/gps/products/%W/igs%W%D.sp3c.Z']

2,'IGS EPH Rapid (CDDIS)',  [cddis,'/gps/products/%W/igr%W%D.sp3.Z']
2,'IGS CLK Rapid (CDDIS)',  [cddis,'/gps/products/%W/igr%W%D.clk.Z']
2,'IGS ERP Rapid (CDDIS)',  [cddis,'/gps/products/%W/igr%W%D.erp.Z']
2,'IGS SUM Rapid (CDDIS)',  [cddis,'/gps/products/%W/igr%W%D.sum.Z']
2,'IGS EPH Rapid(SP3C)(CDDIS)', [cddis,'/gps/products/%W/igr%W%D.sp3c.Z']

2,'IGS EPH URapid (CDDIS)', [cddis,'/gps/products/%W/igu%W%D_%hb.sp3.Z']
2,'IGS ERP URapid (CDDIS)', [cddis,'/gps/products/%W/igu%W%D_%hb.erp.Z']
2,'IGS SUM URapid (CDDIS)', [cddis,'/gps/products/%W/igu%W%D_%hb.sum.Z']
2,'IGS EPH URapid(SP3C)(CDDIS)',[cddis,'/gps/products/%W/igu%W%D_%hb.sp3c.Z']

2,'IGS POS Final (CDDIS)',  [cddis,'/gps/products/%W/igs%yP%W.snx.Z']
2,'IGS POS Final (SOPAC)',  [sopac,'/products/%W/igs%yp%W.snx.Z']
2,'IGS POS Final (IGN)',    [ign,'/igs/products/%W/igs%yP%W.snx.Z']

2,'IGS ZPD Final (CDDIS)',  [cddis,'/gps/products/%W/trop/%s%W.zpd.Z']
2,'IGS ZPD Final (SOPAC)',  [sopac,'/troposphere/%W/%s%W.zpd.Z']
2,'IGS ZPD Final (IGN)',    [ign,'/igs/tropo/%W/%s%W.zpd.Z']
2,'IGS ZPD URapid (CDDIS)', [cddis,'/gps/products/trop/nrt/%W/igs60_%W%D_%g.tro.Z']

2,'IGS TEC Final (CDDIS)',  [cddis,'/gps/products/ionex/%Y/%n/igsg%n0.%yi.Z']
2,'IGS TEC Rapid (CDDIS)',  [cddis,'/gps/products/ionex/%Y/%n/igrg%n0.%yi.Z']
2,'IGS TEC Final (IGN)',    [ign,'/igs/iono/%Y/%n/igsg%n0.%yi.Z']

2,'GSI POS(F1) 電子基準点', [gsi,'/coordinates_F1/%Y/%S.%y.pos']
2,'GSI POS(F2) 電子基準点', [gsi,'/coordinates_F2/%Y/%S.%y.pos']

2,'COD EPH (CDDIS)',        [cddis,'/gps/products/%W/cod%W%D.eph.Z']
2,'COD CLK (CDDIS)',        [cddis,'/gps/products/%W/cod%W%D.clk.Z']
2,'COD ERP (CDDIS)',        [cddis,'/gps/products/%W/cod%W7.erp.Z']
2,'EMR EPH (CDDIS)',        [cddis,'/gps/products/%W/emr%W%D.sp3.Z']
2,'EMR CLK (CDDIS)',        [cddis,'/gps/products/%W/emr%W%D.clk.Z']
2,'EMR ERP (CDDIS)',        [cddis,'/gps/products/%W/emr%W7.erp.Z']
2,'ESA EPH (CDDIS)',        [cddis,'/gps/products/%W/esa%W%D.sp3c.Z']
2,'ESA CLK (CDDIS)',        [cddis,'/gps/products/%W/esa%W%D.clk.Z']
2,'ESA ERP (CDDIS)',        [cddis,'/gps/products/%W/esa%W7.erp.Z']
2,'GFZ EPH (CDDIS)',        [cddis,'/gps/products/%W/gfz%W%D.sp3.Z']
2,'GFZ CLK (CDDIS)',        [cddis,'/gps/products/%W/gfz%W%D.clk.Z']
2,'GFZ ERP (CDDIS)',        [cddis,'/gps/products/%W/gfz%W7.erp.Z']
2,'JPL EPH (CDDIS)',        [cddis,'/gps/products/%W/jpl%W%D.sp3.Z']
2,'JPL CLK (CDDIS)',        [cddis,'/gps/products/%W/jpl%W%D.clk.Z']
2,'JPL ERP (CDDIS)',        [cddis,'/gps/products/%W/jpl%W7.erp.Z']
2,'MIT EPH (CDDIS)',        [cddis,'/gps/products/%W/mit%W%D.sp3.Z']
2,'MIT CLK (CDDIS)',        [cddis,'/gps/products/%W/mit%W%D.clk.Z']
2,'MIT ERP (CDDIS)',        [cddis,'/gps/products/%W/mit%W7.erp.Z']
2,'NGS EPH (CDDIS)',        [cddis,'/gps/products/%W/ngs%W%D.sp3.Z']
2,'NGS ERP (CDDIS)',        [cddis,'/gps/products/%W/ngs%W7.erp.Z']
2,'SIO EPH (CDDIS)',        [cddis,'/gps/products/%W/sio%W%D.sp3.Z']
2,'SIO ERP (CDDIS)',        [cddis,'/gps/products/%W/sio%W7.erp.Z']

2,'COD ION (CODE)',         [code,'/%Y/CODG%n0.%yI.Z']
2,'COD EPH Rapid (CODE)',   [code,'/COD%W%D.EPH_R']
2,'COD CLK Rapid (CODE)',   [code,'/COD%W%D.CLK_R']
2,'COD ERP Rapid (CODE)',   [code,'/COD%W%D.ERP_R']
2,'COD TRO Rapid (CODE)',   [code,'/COD%W%D.TRO_R']
2,'COD ION Rapid (CODE)',   [code,'/CORG%n0.%yI.Z']
2,'COD EPH URapid (CODE)',  [code,'/COD.EPH_U']
2,'COD ERP URapid (CODE)',  [code,'/COD.ERP_U']
2,'COD TRO URapid (CODE)',  [code,'/COD.TRO_U']
2,'COD P1P2 (CODE)',        [code,'/%Y/P1P2%y%m.DCB.Z']
2,'COD P1C1 (CODE)',        [code,'/%Y/P1C1%y%m.DCB.Z']
2,'COD P1P2 ALL (CODE)',    [code,'/P1P2_ALL.DCB']

2,'GRACE LEVEL1B (JPL)',    [podaac,'/grace/data/L1B/%Y/grace_1B_%Y-%m-%d_00.tar.gz']

% eop --------------------------------------------------------------------------

3,'IERS Bulletin A',        'ftp://maia.usno.navy.mil/ser7/ser7.dat'
3,'IERS Bulletin B',        'ftp://hpiers.obspm.fr/iers/bul/bulb/bulletinb.%M'
3,'IERS Bulletin C',        'ftp://hpiers.obspm.fr/iers/bul/bulc/bulletinc.dat'
3,'IERS EOP C04',           'ftp://hpiers.obspm.fr/iers/eop/eopc04/eopc04.%y'

% others -----------------------------------------------------------------------

4,'IGSMAIL',                [cddis,'/gps/igsmail/igsmess.%N']
4,'IGSREPORT',              [cddis,'/gps/igsreport/igsreport.%N']
4,'IGSMAIL(ARCHIVE)',       [cddis,'/gps/igsmail/%Y/igsmess.%N']
4,'IGSREPORT(ARCHIVE)',     [cddis,'/gps/igsreport/%Y/igsreport.%N']
4,'IGSMAIL INDEX',          [cddis,'/gps/igsmail/igsmess_%Y.index']
4,'IGSREPORT INDEX',        [cddis,'/gps/igsreport/igsreport_%Y.index']

% gpv --------------------------------------------------------------------------
4,'JMA MSM SFC T=0-18',     [gfd,'/arch-large2/jmadata/gpv/original/%Y%m%d/MSM%hbSFC018_%{1|2}']
4,'JMA MSM PLM T=0-18',     [gfd,'/arch-large2/jmadata/gpv/original/%Y%m%d/MSM%hbPLM018_%{1|2}']
4,'JMA MSM PMH T=0-18',     [gfd,'/arch-large2/jmadata/gpv/original/%Y%m%d/MSM%hbPMH018']
4,'JMA RSM SFC T=0-24',     [gfd,'/arch-large2/jmadata/gpv/original/%Y%m%d/RSM%hcSFC024']
4,'JMA RSM SFC T=25-51',    [gfd,'/arch-large2/jmadata/gpv/original/%Y%m%d/RSM%hcSFC051']
4,'JMA RSM PLM T=0-24',     [gfd,'/arch-large2/jmadata/gpv/original/%Y%m%d/RSM%hcPLM024']
4,'JMA RSM PLM T=27-51',    [gfd,'/arch-large2/jmadata/gpv/original/%Y%m%d/RSM%hcPLM051']
4,'JMA RSM PMH T=0-24',     [gfd,'/arch-large2/jmadata/gpv/original/%Y%m%d/RSM%hcPMH024']
4,'JMA RSM PMH T=27-51',    [gfd,'/arch-large2/jmadata/gpv/original/%Y%m%d/RSM%hcPMH051']
4,'JMA GSM T=0-24',         [gfd,'/arch-large2/jmadata/gpv/original/%Y%m%d/GSM%hcX024']
4,'JMA GSM T=30-48',        [gfd,'/arch-large2/jmadata/gpv/original/%Y%m%d/GSM%hcX048']
4,'JMA GSM T=54-84',        [gfd,'/arch-large2/jmadata/gpv/original/%Y%m%d/GSM%hcX084']
4,'JMA GSM T=96-192',       [gfd,'/arch-large2/jmadata/gpv/original/%Y%m%d/GSM%hcX192']
4,'JMA GSM T=90-180',       [gfd,'/arch-large2/jmadata/gpv/original/%Y%m%d/GSM%hcX180']
4,'JMA GSM T=90-180',       [gfd,'/arch-large2/jmadata/gpv/original/%Y%m%d/GSM%hcX180']
4,'JMA GSM N PRS T=0-72(DDB)',[jma,'/DATA/jp034/g003f12/%y%m%d12/hpx%{a|e|i|k}89']
4,'JMA GSM N HEI T=0-72(DDB)',[jma,'/DATA/jp034/g003f12/%y%m%d12/hhx%{a|e|i|k}%{85|70|50|30|25|20|10}']
4,'JMA GSM N TEM T=0-72(DDB)',[jma,'/DATA/jp034/g003f12/%y%m%d12/htx%{a|e|i|k}%{98|85|70|50|30|25|20|10}']
4,'JMA GSM N WIN T=0-72(DDB)',[jma,'/DATA/jp034/g003f12/%y%m%d12/h%{u|v}x%{a|e|i|k}%{98|85|70|50|30|25|20|10}']
4,'JMA GSM N DEW T=0-72(DDB)',[jma,'/DATA/jp034/g003f12/%y%m%d12/hrx%{a|e|i|k}%{85|70}']
4,'JMA GSM S PRS T=0-72(DDB)',[jma,'/DATA/jp034/g004f12/%y%m%d12/hpy%{a|e|i|k}89']
4,'JMA GSM S HEI T=0-72(DDB)',[jma,'/DATA/jp034/g004f12/%y%m%d12/hhy%{a|e|i|k}%{85|70|50|30|25|20|10}']
4,'JMA GSM S TEM T=0-72(DDB)',[jma,'/DATA/jp034/g004f12/%y%m%d12/hty%{a|e|i|k}%{98|85|70|50|30|25|20|10}']
4,'JMA GSM S WIN T=0-72(DDB)',[jma,'/DATA/jp034/g004f12/%y%m%d12/h%{u|v}y%{a|e|i|k}%{98|85|70|50|30|25|20|10}']
4,'JMA GSM S DEW T=0-72(DDB)',[jma,'/DATA/jp034/g004f12/%y%m%d12/hry%{a|e|i|k}%{85|70}']
};
