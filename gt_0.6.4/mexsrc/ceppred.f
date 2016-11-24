CTITLE CEPPRED
*     This program has been modifed slightly from Herring's KSV
*     program to account for the bias in celestial pole offsets
*     (-43.1 in dpsi and -5.1 in depsilon).  These bias offsets
*     are added to the results in the calculation of dpsi_tot
*     and deps_tot.  These biases are chosen to best fit the
*     model to the VLBI observations.
      subroutine KSV_1996_3( jd, dpsi_ls, deps_ls,
     .                    dpsi_plan, deps_plan,
     .                    dpsi_fcn , deps_fcn ,
     .                    dpsi_prec, deps_prec,
     .                    dpsi_tot , deps_tot    )
 
*     Subroutine to compute the complete KSV_1996_3 nutation series
*     with the associated corrections for planetary nutations,
*     the freely excited free-core-nutation (valid for 1988-1995),
*     the precession constant change and a rate of change of oblquity.)
 
* USAGE:
*     call  KSV_1996_3( jd, dpsi_ls, deps_ls, dpsi_plan, deps_plan,
*    .        dpsi_fcn , deps_fcn , dPsi_prec, deps_prec,
*    .       dpsi_tot , deps_tot )
 
* where:
*     <jd>    is the full julian date including fractional part of
*             of the day (REAL*8 input)
*     <dpsi_ls> and <deps_ls> are the luni-solar nutation in
*             longitude and oblquity (mas) (REAL* OUTPUT)
*     <dpsi_plan> and <deps_plan> are the contributions to the
*             nutations in longitude and obliquity due direct
*             planetary nutations and the perturbations of the
*             lunar and terrestrial orbits (mas). (REAL* OUTPUT)
*     <dpsi_fcn> and <deps_fcn> are the contributions to the
*             nutations in longitude and obliquity due the free-
*             excitation of the Free-core-nutation (mas).  These
*             values are valid for 1988-1995. (REAL* OUTPUT)
*     <dpsi_prec> and <deps_prec> are the contributions to the
*             nutations in longitude and obliquity due changes in
*             the precession constant and rate of change of
*             obliquity (mas) (REAL* OUTPUT).
*     <dpsi_tot> and <deps_tot> are the total nutations in longitude
*             and obliquity including the correction for the precession
*             constant (when precession is computed using the IAU 1976
*             precession constant), and are obtained by summing all
*             of the above corrections (mas) (REAL* OUTPUT).
 
* RESTRICTIONS: if <jd> is less than 2000000.0 this routine
*               assumes an MJD has been passed and the time
*               used will be converted to JD.  A warning
*               message will be printed.  See individual modules
*               for further restrictions.
 
* PASSED VARIABLES
*
* INPUT Values
* jd     - Time at which value needed. (jd + fraction of day)
 
* OUTPUT Values
*     dpsi_ls and deps_ls      - luni-solar nutation in
*             longitude and oblquity (mas) (REAL* OUTPUT)
*     dpsi_plan and deps_plan  - contributions to the
*             nutations in longitude and obliquity due direct
*             planetary nutations and the perturbations of the
*             lunar and terrestrial orbits (mas). (REAL* OUTPUT)
*     dpsi_fcn and deps_fcn    - contributions to the
*             nutations in longitude and obliquity due the free-
*             excitation of the Free-core-nutation (mas).  These
*             values are valid for 1988-1994. (REAL* OUTPUT)
*     dpsi_prec and deps_prec  - contributions to the
*             nutations in longitude and obliquity due changes in
*             the precession constant and rate of change of
*             obliquity (mas) (REAL* OUTPUT).
*     dpsi_tot and deps_tot    - total nutations in longitude
*             and obliquity including the correction for the precession
*             constant (when precession is computed using the IAU 1976
*             precession constant), and are obtained by summing all
*             of the above corrections (mas) (REAL* OUTPUT).
 
 
      real*8 jd, dpsi_ls, deps_ls, dpsi_plan, deps_plan,
     .    dpsi_fcn ,  deps_fcn, dpsi_prec, deps_prec,
     .    dpsi_tot, deps_tot
 
*---------------------------------------------------------------
 
*     Call each of the routines needed for each contribution.
 
*     Luni-solar nutation
      call ls_nut( jd, dpsi_ls, deps_ls )
 
*     Planetary nutation
      call plan_nut ( jd, dpsi_plan, deps_plan )
 
*     Freely excited FCN (NOTE: No warning message is printed
*     if the JD is out of the range of 1988-1994)
      call fcn_nut ( jd, dpsi_fcn , deps_fcn )
 
*     Precession and obliquity rate contributions (NOTE: IAU-1976
*     precession constant assumed to be used in the basic calculation
*     of precession).
 
      call prec_nut( jd, dpsi_prec, deps_prec )
 
*     Now add up all of the terms to get the total nutation angles
*     and add on the bias offsets.
 
      dpsi_tot = dpsi_ls + dpsi_plan + dpsi_fcn + dpsi_prec - 43.1d0
      deps_tot = deps_ls + deps_plan + deps_fcn + deps_prec - 5.1d0
 
*     Thats all
      return
      end
 
CTITLE LS_NUT
 
      subroutine ls_nut( jd, dpsi_ls, deps_ls )
 
*     Routine to compute the KSV_1996_3 luni-solar contributions
*     the nutations in longitude and obliquity.  The KSV_1996_3 is
*     based on:
 
*     (1) The Souchay and Kinoshita Rigid Earth nutation series
*     KSRE95.   The seven terms with duplicate
*     arguments on the KS series have been combined into single terms.
*     The arguments were:
*      l   lp  F   D  Om  Period (days)
*     -1   0   0   1   0   411.78
*      0   1   0   0   0   365.26
*      0  -1   2  -2   2   365.22
*      0   2   0   0   0   182.63
*      0   0   2  -2   2   182.62
*      0  -1  -2   2  -2   121.75
*      0   0   0   1   0    29.53
*     The series has also been sorted in increasing order of amplitude
*     of the nutation.  This is to minimize rounding error as the series
*     is summed in reverse order.
 
*     (2) Estimates of the Retrograde FCN resonance factors from
*     the Mathews et al., nutation formulation (full complex
*     estimates, and scaling parameter R from the same
*     theory.  
*     (The definition on complex values here is Rr + iRi; in earlier
*      versions the definition was Rr -iRi.  This change does not effect
*      the nutation series coefficients).  The nutation amplitudes are
*      still defined as ar - ai).
*
*   The resonance factors used are:
*   Type             Amplitude                    Frequency (cpsd)
*                Real            Imag.         Real            Imag.
*   RFCN     -.00011489752   .00000214130     -1.00231888314  -.00002920327 
*   CW       -.00057992000   .00000000000       .00253170000   .00000000000  
*   R and R' 1.04901828112  -.00150732471      -.25517427386   .03965769073

*     (3) The effects of annual modulation of geodetic precession.
*     The correction applied is
*      0  1  0  0  0 -0.150 (correction to in-phase nutation in
*                         longitude).

*     (4) A prograde annual nutation has been estimated along with the
*      resonance coefficients.  This probably reflects the influence of
*      S1 atmospheric tide.

*     (5) The free RFCN mode was estimated once every two years for the
*      data after 1984.  (See values commented in eval_ls_nut.  For the
*      last 6 years the values seem to be resonably stable.  

*     (6) The new Simons et al., fundamental arguments are used in this
*      version.  (The largest change from KSV_1995_1 was 0.007 mas for
*      semiannual nutation.  All other changes, including the 18.6 year
*      nutation were 0.001-0.002 mas.)

 
*     REFERENCES:
* NEW Version based on: Corrections and new developments in rigid Earth
*     nutation theory: Lunisolar influence including indirect planetary
*     effects, J. Souchay and H. Kinioshita, Astron. and Astrophys., 1995.
* (Version here based on data files: KSRE95_FIG_PSI.DAT KSRE95_FIG_EPS.DAT
*  and generated with ks_plan.f)
*     Kinoshita, H., and J. Souchay, The theory of the nutations for
*         the rigid Earth at the second order, Celes. Mech. and Dynam.
*         Astron., 48, 187--266, 1990.
*     Mathews, P. M., B. A. Buffett, T. A. Herring, and I. I. Shapiro,
*         Forced nutations of the Earth: Influence of the inner core
*         dynamics, 1, Theory. J. Geophs. Res., 96, 8219--8242, 1991.
*     Simon, J. L., Bretagnon, P., Chapront, J., Chapront-Touze, M.,
*         Francou, G., Laskar, J., 1994, "Numerical Expressions for 
*         Precession Formulae and Mean Elements for the Moon and
*         Planets," Astron. Astrophys., 282, pp. 663-683.

 
 
* USAGE:
*     call ls_nut( jd, dpsi_ls, deps_ls )
*     where <jd>    is a full julian date with fractional part
*                   of the day added (REAL*8 INPUT)
*     and <dpsi_ls> and <deps_ls> are the nutations
*                   in longitude and obliquity in milliarcsec.
*                   (REAL*8 OUTPUT)
 
* RESTRICTIONS: if <jd> is less than 2000000.0 this routine
*               assumes an MJD has been passed and the time
*               used will be converted to JD.  A warning
*               message will be printed.
 
* PASSED VARIABLES
*
* INPUT Values
* jd     - Time at which value needed. (jd + fraction of day)
 
* OUTPUT Values
* dpsi_ls  - The nutation in longitude (mas).
* deps_ls  - The nutation in obliquity (mas).
 
 
 
      real*8 jd, dpsi_ls, deps_ls
 
* LOCAL VARIABLES
 
*   epoch       - Julian date (jd passed in unless the JD
*                 appears to be an MJD in which case it is
*                 converted to JD (2 400 000.5d0 added)
*   ls_arg(5)   - The arguments for the Luni-solar nutations.
*                 (l, l', F, D and Omega).  All in Radians.
 
 
      real*8 epoch, ls_arg(5)
 
***** Check to make sure user passed JD and not MJD.  Correct
*     problem and warn the user.
      if( jd .lt.2 000 000.0d0  ) then
          write(*,100) jd
 100      format('**WARNING** MJD apparently passed to SD_COMP',
     .          ' Value (',F10.2,') converted to JD')
          epoch = jd + 2 400 000.5d0
      else
          epoch = jd
      end if
 
***** Get the fundamental arguments at this epoch
 
      call ls_angles( epoch, ls_arg)
 
*     Now compute the luni-solare nutations by summing over all
*     terms in the series.
 
      call eval_ls_nut( epoch, ls_arg, dpsi_ls, deps_ls )
 
*     Thats all
      return
      end
 
 
CTITLE LS_ANGLES
 
      subroutine ls_angles( epoch, ls_arg )
 
*     Routine to compute the value of the fundamental argument
*     for Brown's arguments.  Arguments based on the IERS
*     standards.

* MOD TAH 960206: Changed arguments to use Simons et al., 1994 
* values:
*     Simon, J. L., Bretagnon, P., Chapront, J., Chapront-Touze, M.,
*          Francou, G., Laskar, J., 1994, "Numerical Expressions for 
*          Precession Formulae and Mean Elements for the Moon and
*          Planets," Astron. Astrophys., 282, pp. 663-683.


* PHYSICAL CONSTANTS
 
*   pi          - Define here to full precision
*   rad_to_deg  - Conversion from radians to degs.
*   DJ2000      - Julian date of J2000
*   sec360      - number of seconds in 360 degreees.
 
 
      real*8 pi, rad_to_deg, DJ2000, sec360
 
      parameter ( pi            = 3.1415926535897932D0 )
      parameter ( DJ2000        = 2451545.d0           )
      parameter ( sec360        = 1296000.d0           )
 
*     Computed quanities
      parameter ( rad_to_deg    = 180.d0   /pi         )
 
*-------------------------------------------------------------------
 
* PASSED VARIABLES
 
* INPUT
* epoch  - Julian date for arguments (jd + fraction of day, REAL*8)
 
* OUTPUT
* ls_arg(5) -  Brown's arguments (radians, REAL*8)
 
 
      real*8 epoch, ls_arg(5)
 
* LOCAL VARIABLES
*      cent             - Julian centuries to DJ2000.
*      el,eld           - Mean longitude of moon minus mean
*                       - longitude of moon's perigee (arcsec)
*      elc(5)           - Coefficients for computing el
*      elp,elpd         - Mean longitude of the sun minus mean
*                       - longitude of sun perigee (arcsec)
*      elpc(5)          - Coeffiecents for computing elp
*      f,fd             - Moon's mean longitude minus omega (sec)
*      fc(5)            - Coefficients for computing f
*      d,dd             - Mean elongation of the moon from the
*                       - sun (arcsec)
*      dc(5)            - coefficients for computing d
*      om,omd           - longitude of the ascending node of the
*                       - moon's mean orbit on the elliptic
*                       - measured from the mean equinox of date
*      omc(5)           - Coefficients for computing om.
 
 
      real*8 cent, el,eld, elc(5), elp, elpd, elpc(5),
     .    f,fd, fc(5), d,dd, dc(5), om,omd, omc(5)
 
****  DATA statements for the fundamental arguments.
*     Simons et al., 1994 values
 
      data elc    /    -0.00024470d0,    0.051635d0,   31.8792d0,    
     .         1717915923.2178d0,   485868.249036d0/
      data elpc   /    -0.00001149d0,    -0.000136d0,  -0.5532d0,   
     .          129596581.0481d0,   1287104.79305d0/
      data fc     /     0.00000417d0,    -0.001037d0,  -12.7512d0,
     .         1739527262.8478d0,    335779.526232d0/
      data dc     /    -0.00003169d0,     0.006593d0,   -6.3706d0,   
     .         1602961601.2090d0,   1072260.70369d0/
* MOD TAH KSV_1996_3: 960606: Replaced <Om> with expression from b.3 of 
*     Simon et al., 1994 since b.3 is computed with new precession constant
*     (Only the rate changes).   
      data omc    /    -0.00005939,       0.007702d0,    7.4722d0, 
     .           -6962890.5431d0,     450160.398036d0/

 
****  Get the number of centuries to current time
 
      cent = (epoch-dj2000) / 36525.d0
 
****  Compute angular arguments and their time derivatives
* New formulas adding in the higher order term.

      el = elc(1) * cent**4 + elc(2) * cent**3 + elc(3) * cent**2
     .      + elc(4) * cent + elc(5)
      el = mod( el, sec360 )
      eld = 4.d0 * elc(1) * cent**3 + 3.d0 * elc(2) * cent**2 + 
     .      2.d0 * elc(3) * cent    +        elc(4) 
c
      elp = elpc(1) * cent**4 + elpc(2) * cent**3 + elpc(3) * cent**2
     .     + elpc(4) * cent + elpc(5)
      elp = mod( elp, sec360 )
      elpd = 4.d0 * elpc(1) * cent**3 + 3.d0 * elpc(2) * cent**2 + 
     .       2.d0 * elpc(3) * cent    +        elpc(4) 
c
      f = fc(1) * cent**4 + fc(2) * cent**3 + fc(3) * cent**2
     .     + fc(4) * cent + fc(5)
      f = mod( f, sec360 )
      fd = 4.d0 * fc(1) * cent**3 + 3.d0 * fc(2) * cent**2 + 
     .     2.d0 * fc(3) * cent    +        fc(4) 
c
      d = dc(1) * cent**4 + dc(2) * cent**3 + dc(3) * cent**2
     .   + dc(4) * cent + dc(5)
      d = mod( d, sec360 )
      dd = 4.d0 * dc(1) * cent**3 + 3.d0 * dc(2) * cent**2 + 
     .     2.d0 * dc(3) * cent    +        dc(4) 
c
      om = omc(1) * cent**4 + omc(2) * cent**3 + omc(3) * cent**2
     .     + omc(4) * cent + omc(5)
      om = mod( om, sec360 )
      omd = 4.d0 * omc(1) * cent**3 + 3.d0 * omc(2) * cent**2 + 
     .      2.d0 * omc(3) * cent    +        omc(4) 
c
c
 
****  Now save the values.  Convert values from arcseconds to radians
 
      ls_arg(1) = el / (3600.d0*rad_to_deg)
      ls_arg(2) = elp/ (3600.d0*rad_to_deg)
      ls_arg(3) = f  / (3600.d0*rad_to_deg)
      ls_arg(4) = d  / (3600.d0*rad_to_deg)
      ls_arg(5) = om / (3600.d0*rad_to_deg)
 
***** Thats all
      return
      end
 
CTITLE EVAL_LS_NUT
 
      subroutine eval_ls_nut( epoch, ls_arg, dpsi_ls, deps_ls )
 
*     Routine to compute the nutations in longitude and obliquity
*     by summing over all terms in the nutations series.
 
* NOTE: ls_angles must be called before routine.
 
* PARAMETERS:
 
* num_ls  - Number of terms in the nutations series
 
      integer*4 num_ls
 
      parameter ( num_ls      =  263)
 
*   DJ2000      - Julian date of J2000
*   pi          - Pi.
 
 
      real*8 pi, DJ2000
 
      parameter ( pi            = 3.1415926535897932D0 )
      parameter ( DJ2000        = 2451545.d0           )
 
* PASSED PARAMETERS:
 
* INPUT:
* epoch      - Julian date at which nutation angles are needed.
* ls_arg(5)  - Five arguments for the nutatoions (l,l',F,D and Om)
*              computed at the epoch that the nutations need to be
*              evaluated (rad) (REAL*8)
 
* OUTPUT:
* dpsi_ls, deps_ls   - nutations in longitude and obliquity (mas)
*              (REAL*8)
 
 
      real*8 epoch, ls_arg(5), dpsi_ls, deps_ls
 
* LOCAL VARIABLES:
 
*  i and j   - Counters to loop over the coeffients and the argumemts
 
 
      integer*4 i,j
 
*  arg       - Final summed argumemt for the nutations
*              contributions (rads)
*  cent      - Number of centuries since J2000.
*  dpsi_lsu and deps_lsu - Nutations in longitude and oblquity
*              in micro-arc-sec (units that the data statements
*              are in)
 
 
      real*8 arg, cent, dpsi_lsu, deps_lsu
 
*     Fortan data statements for the nutation series
*     Series based on the following model:
*     Type                Amplitude           Frequency
*       1     -.00011489814   .00000214165     -1.00231887299  -.00002920212
*       2      .00000000000   .00000000000      -.99780880000   .00000000000
*       3     -.00057992000   .00000000000       .00253170000   .00000000000
*       4      .00000000000   .00000000000       .00041390000   .00000000000
*       5     1.04901875266  -.00150730877      -.25518016085   .03966200327

*       6      .00000000000   .00000000000       .00000000000   .00000000000
*    
*     RFCN Freq.  -1.00231887299 cyc per sidreal day, Period   430.0665 solar days
*  IX01-IX27(11,10)  -- Invidual declarations of the coefficents
*                           of the nutation series so no data statement has more than 10 lines
*                           The first 5 values are the arguments for l lp F D Om
*                           The remaining elements are:
*                            6 - Nutation in longitude psi (sin, uas)
*                            7 - dpsi/dt (uasec/cent)
*                            8 - Nutation in oblquity eps (cos, uas)
*                            9 - deps/dt (uas/cent)
*                           10 - Out-of-phase longitude (cos, uas)
*                           11 - Out-of-phase obliquity (sin, uas)

      integer*4 IX01(11,10), IX02(11,10), IX03(11,10), IX04(11,10), 
     .          IX05(11,10), IX06(11,10), IX07(11,10), IX08(11,10), 
     .          IX09(11,10), IX10(11,10), IX11(11,10), IX12(11,10), 
     .          IX13(11,10), IX14(11,10), IX15(11,10), IX16(11,10), 
     .          IX17(11,10), IX18(11,10), IX19(11,10), IX20(11,10), 
     .          IX21(11,10), IX22(11,10), IX23(11,10), IX24(11,10), 
     .          IX25(11,10), IX26(11,10), IX27(11, 3)

      integer*4 nutc_int(11,263)
      equivalence (nutc_int(1,  1),IX01(1,1))
      equivalence (nutc_int(1, 11),IX02(1,1))
      equivalence (nutc_int(1, 21),IX03(1,1))
      equivalence (nutc_int(1, 31),IX04(1,1))
      equivalence (nutc_int(1, 41),IX05(1,1))
      equivalence (nutc_int(1, 51),IX06(1,1))
      equivalence (nutc_int(1, 61),IX07(1,1))
      equivalence (nutc_int(1, 71),IX08(1,1))
      equivalence (nutc_int(1, 81),IX09(1,1))
      equivalence (nutc_int(1, 91),IX10(1,1))
      equivalence (nutc_int(1,101),IX11(1,1))
      equivalence (nutc_int(1,111),IX12(1,1))
      equivalence (nutc_int(1,121),IX13(1,1))
      equivalence (nutc_int(1,131),IX14(1,1))
      equivalence (nutc_int(1,141),IX15(1,1))
      equivalence (nutc_int(1,151),IX16(1,1))
      equivalence (nutc_int(1,161),IX17(1,1))
      equivalence (nutc_int(1,171),IX18(1,1))
      equivalence (nutc_int(1,181),IX19(1,1))
      equivalence (nutc_int(1,191),IX20(1,1))
      equivalence (nutc_int(1,201),IX21(1,1))
      equivalence (nutc_int(1,211),IX22(1,1))
      equivalence (nutc_int(1,221),IX23(1,1))
      equivalence (nutc_int(1,231),IX24(1,1))
      equivalence (nutc_int(1,241),IX25(1,1))
      equivalence (nutc_int(1,251),IX26(1,1))
      equivalence (nutc_int(1,261),IX27(1,1))
      data IX01/   0,   0,   0,   0,   1, -17206262, -17419,
     .                    9205348,  886, 3608, 1543,
     .            0,   0,   2,  -2,   2,  -1317015,   -156,
     .                     573059, -306,-1398, -463,
     .            0,   0,   2,   0,   2,   -227718,    -23,
     .                      97864,  -48,  288,  145,
     .            0,   0,   0,   0,   2,    207428,     21,
     .                     -89746,   47,  -71,  -29,
     .            0,   1,   0,   0,   0,    147545,   -364,
     .                       7390,  -19, 1121, -198,
     .            0,   1,   2,  -2,   2,    -51687,    123,
     .                      22440,  -68,  -54,  -18,
     .            1,   0,   0,   0,   0,     71118,      7,
     .                       -687,    0,  -98,   39,
     .            0,   0,   2,   0,   1,    -38752,    -37,
     .                      20076,    2,   37,   34,
     .            1,   0,   2,   0,   2,    -30136,     -4,
     .                      12896,   -6,   81,   37,
     .            0,  -1,   2,  -2,   2,     21583,    -49,
     .                      -9591,   30,    6,   12 /
      data IX02/   0,   0,   2,  -2,   1,     12820,     14,
     .                      -6897,   -1,   18,    4,
     .           -1,   0,   2,   0,   2,     12353,      1,
     .                      -5333,    3,    1,   -1,
     .           -1,   0,   0,   2,   0,     15700,      1,
     .                       -127,    0,  -19,    9,
     .            1,   0,   0,   0,   1,      6314,      6,
     .                      -3323,    0,    3,   -1,
     .           -1,   0,   0,   0,   1,     -5797,     -6,
     .                       3141,    0,  -19,   -8,
     .           -1,   0,   2,   2,   2,     -5965,     -1,
     .                       2554,   -1,   15,    7,
     .            1,   0,   2,   0,   1,     -5163,     -4,
     .                       2635,    0,   12,    8,
     .           -2,   0,   2,   0,   1,      4590,      5,
     .                      -2424,   -1,    1,    1,
     .            0,   0,   0,   2,   0,      6336,      1,
     .                       -125,    0,  -16,    3,
     .            0,   0,   2,   2,   2,     -3854,      0,
     .                       1643,    0,   15,    7 /
      data IX03/  -2,   0,   0,   2,   0,     -4774,      0,
     .                         48,    0,   -2,   -3,
     .            2,   0,   2,   0,   2,     -3102,      0,
     .                       1322,   -1,   13,    6,
     .            1,   0,   2,  -2,   2,      2863,      0,
     .                      -1234,    1,    0,    0,
     .           -1,   0,   2,   0,   1,      2044,      2,
     .                      -1076,    0,    1,    0,
     .            2,   0,   0,   0,   0,      2923,      0,
     .                        -62,    0,   -8,    1,
     .            0,   0,   2,   0,   0,      2585,      0,
     .                        -56,    0,   -7,    1,
     .            0,   1,   0,   0,   1,     -1406,     -3,
     .                        857,    0,    8,   -4,
     .           -1,   0,   0,   2,   1,      1517,      1,
     .                       -801,    0,    1,    0,
     .            0,   2,   2,  -2,   2,     -1578,      7,
     .                        685,   -4,   -2,   -1,
     .            0,   0,  -2,   2,   0,      2178,      0,
     .                        -15,    0,    1,    1 /
      data IX04/   1,   0,   0,  -2,   1,     -1286,     -1,
     .                        694,    0,   -4,   -2,
     .            0,  -1,   0,   0,   1,     -1269,      1,
     .                        642,    1,    6,    2,
     .           -1,   0,   2,   2,   1,     -1022,     -1,
     .                        522,    0,    2,    2,
     .            0,   2,   0,   0,   0,      1671,     -8,
     .                         14,    0,   -1,    1,
     .            1,   0,   2,   2,   2,      -768,      0,
     .                        325,    0,    4,    2,
     .           -2,   0,   2,   0,   0,     -1102,      0,
     .                         10,    0,   -1,    0,
     .            0,   1,   2,   0,   2,       757,     -2,
     .                       -326,   -2,   -1,   -1,
     .            0,   0,   2,   2,   1,      -664,     -1,
     .                        335,   -1,    2,    1,
     .            0,  -1,   2,   0,   2,      -714,      2,
     .                        307,    2,    1,    0,
     .            0,   0,   0,   2,   1,      -631,     -1,
     .                        327,    0,    0,    0 /
      data IX05/   1,   0,   2,  -2,   1,       580,      1,
     .                       -307,    0,    0,    0,
     .            2,   0,   2,  -2,   2,       643,      0,
     .                       -277,    0,   -1,    0,
     .           -2,   0,   0,   2,   1,      -579,     -1,
     .                        304,    0,   -1,    0,
     .            2,   0,   2,   0,   1,      -533,      0,
     .                        269,    0,    2,    1,
     .            0,  -1,   2,  -2,   1,      -477,     -1,
     .                        271,   -1,    0,    0,
     .            0,   0,   0,  -2,   1,      -493,     -1,
     .                        272,    0,   -2,   -1,
     .           -1,  -1,   0,   2,   0,       735,      0,
     .                         -5,    0,   -1,    0,
     .            2,   0,   0,  -2,   1,       405,      0,
     .                       -220,    0,    1,    0,
     .            1,   0,   0,   2,   0,       657,      0,
     .                        -20,    0,   -2,    0,
     .            0,   1,   2,  -2,   1,       361,      0,
     .                       -194,    0,    1,    0 /
      data IX06/   1,  -1,   0,   0,   0,       471,      0,
     .                         -4,    0,   -1,    0,
     .           -2,   0,   2,   0,   2,      -311,      0,
     .                        131,    0,    0,    0,
     .            3,   0,   2,   0,   2,      -289,      0,
     .                        124,    0,    2,    1,
     .            0,  -1,   0,   2,   0,       435,      0,
     .                         -9,    0,   -1,    0,
     .            1,  -1,   2,   0,   2,      -287,      0,
     .                        123,    0,    1,    0,
     .           -1,  -1,   2,   2,   2,      -282,      0,
     .                        122,    0,    1,    0,
     .            0,   0,   0,   1,   0,      -422,      0,
     .                          3,    0,    1,    0,
     .           -1,   0,   2,   0,   0,      -404,      0,
     .                          4,    0,    1,    0,
     .            0,  -1,   2,   2,   2,      -264,      0,
     .                        114,    0,    1,    0,
     .           -2,   0,   0,   0,   1,      -228,      0,
     .                        126,    0,   -1,    0 /
      data IX07/   1,   1,   2,   0,   2,       246,      0,
     .                       -106,    0,   -1,    0,
     .            2,   0,   0,   0,   1,       218,      0,
     .                       -114,    0,    0,    0,
     .           -1,   1,   0,   1,   0,       327,      0,
     .                         -1,    0,    0,    0,
     .            1,   1,   0,   0,   0,      -338,      0,
     .                          4,    0,    0,    0,
     .            1,   0,   2,   0,   0,       334,      0,
     .                        -11,    0,   -1,    0,
     .           -1,   0,   2,  -2,   1,      -199,      0,
     .                        107,    0,   -1,    0,
     .            1,   0,   0,   0,   2,      -197,      0,
     .                         85,    0,    0,    0,
     .           -1,   0,   0,   1,   0,       405,      0,
     .                        -55,    0,  -35,  -14,
     .            0,   0,   2,   1,   2,       165,      0,
     .                        -72,    0,    0,    0,
     .           -1,   0,   2,   4,   2,      -151,      0,
     .                         66,    0,    1,    0 /
      data IX08/   0,  -2,   2,  -2,   1,      -130,      0,
     .                         69,    0,    0,    0,
     .           -1,   1,   0,   1,   1,       132,      0,
     .                        -68,    0,    0,    0,
     .            1,   0,   2,   2,   1,      -133,      0,
     .                         66,    0,    1,    0,
     .           -2,   0,   2,   2,   2,       139,      0,
     .                        -60,    0,    0,    0,
     .           -1,   0,   0,   0,   2,       139,      0,
     .                        -60,    0,    0,    0,
     .            1,   1,   2,  -2,   2,       128,      0,
     .                        -55,    0,    0,    0,
     .           -2,   0,   2,   4,   2,      -121,      0,
     .                         52,    0,    0,    0,
     .           -1,   0,   4,   0,   2,       115,      0,
     .                        -49,    0,    0,    0,
     .            2,   0,   2,  -2,   1,       101,      0,
     .                        -54,    0,    0,    0,
     .            2,   0,   2,   2,   2,      -108,      0,
     .                         47,    0,    1,    0 /
      data IX09/   1,   0,   0,   2,   1,       -95,      0,
     .                         49,    0,    0,    0,
     .            3,   0,   0,   0,   0,       157,      0,
     .                         -5,    0,   -1,    0,
     .            3,   0,   2,  -2,   2,        94,      0,
     .                        -40,    0,    0,    0,
     .            0,   0,   4,  -2,   2,        91,      0,
     .                        -39,    0,    0,    0,
     .            0,   0,  -2,   2,   1,        87,      0,
     .                        -44,    0,    0,    0,
     .            0,   1,   2,   0,   1,        81,      0,
     .                        -42,    0,    0,    0,
     .            0,   0,   2,  -2,   3,       123,      0,
     .                        -20,    0,    0,    0,
     .           -1,   0,   0,   4,   0,       133,      0,
     .                         -4,    0,    0,    0,
     .            2,   0,  -2,   0,   1,        71,      0,
     .                        -38,    0,    0,    0,
     .           -2,   0,   0,   4,   0,       128,      0,
     .                          1,    0,    0,    0 /
      data IX10/  -1,  -1,   0,   2,   1,        75,      0,
     .                        -39,    0,    0,    0,
     .           -2,  -1,   0,   2,   0,      -115,      0,
     .                          1,    0,    0,    0,
     .            0,  -1,   2,   0,   1,       -66,      0,
     .                         35,    0,    0,    0,
     .           -1,   0,   0,   1,   1,       101,      0,
     .                        -49,    0,   -3,   -1,
     .            0,   0,  -2,   0,   1,       -68,      0,
     .                         36,    0,    0,    0,
     .            0,   1,   0,   0,   2,        69,      0,
     .                        -33,    0,   -1,    0,
     .            0,   0,   2,  -1,   2,       -74,      0,
     .                         31,    0,    0,    0,
     .            0,   0,   2,   4,   2,       -69,      0,
     .                         29,    0,    0,    0,
     .            1,   1,   0,  -2,   1,       -61,      0,
     .                         32,    0,    0,    0,
     .           -1,   1,   0,   2,   0,       -94,      0,
     .                          0,    0,    0,    0 /
      data IX11/   1,  -1,   2,   2,   2,       -59,      0,
     .                         25,    0,    0,    0,
     .            1,  -1,   0,   0,   1,        51,      0,
     .                        -27,    0,    0,    0,
     .            0,   1,  -2,   2,   0,       -90,      0,
     .                          3,    0,    0,    0,
     .            3,   0,   2,   0,   1,       -50,      0,
     .                         25,    0,    0,    0,
     .           -1,   1,   2,   2,   2,        56,      0,
     .                        -24,    0,    0,    0,
     .            0,   1,   2,   2,   2,        54,      0,
     .                        -22,    0,    0,    0,
     .           -1,   0,   0,  -2,   1,       -50,      0,
     .                         27,    0,    0,    0,
     .           -1,   1,   0,   1,   2,       -52,      0,
     .                         23,    0,    0,    0,
     .            0,  -1,   2,   2,   1,       -44,      0,
     .                         24,    0,    0,    0,
     .            1,   0,   2,  -4,   1,       -47,      0,
     .                         24,    0,    0,    0 /
      data IX12/  -1,   0,  -2,   2,   0,        77,      0,
     .                          0,    0,    0,    0,
     .           -1,  -1,   2,   2,   1,       -46,      0,
     .                         24,    0,    0,    0,
     .            0,  -1,   0,   0,   2,        59,      0,
     .                        -25,    0,    0,    0,
     .            2,  -1,   2,   0,   2,       -48,      0,
     .                         21,    0,    0,    0,
     .            1,  -1,   2,   0,   1,       -42,      0,
     .                         22,    0,    0,    0,
     .            0,   0,   0,   2,   2,       -46,      0,
     .                         20,    0,    0,    0,
     .            0,   1,   0,   2,   0,       -67,      0,
     .                          0,    0,    0,    0,
     .           -1,   1,   2,   0,   2,        47,      0,
     .                        -20,    0,    0,    0,
     .            0,   3,   2,  -2,   2,       -44,      0,
     .                         19,    0,    0,    0,
     .            0,  -1,  -2,   2,   0,        66,      0,
     .                          0,    0,    0,    0 /
      data IX13/   0,   0,   0,   1,   1,       -37,      0,
     .                         20,    0,    0,    0,
     .           -1,   0,   2,   2,   0,        64,      0,
     .                          1,    0,    0,    0,
     .            1,   1,   2,   0,   1,        36,      0,
     .                        -18,    0,    0,    0,
     .            2,   1,   2,   0,   2,        40,      0,
     .                        -17,    0,    0,    0,
     .            0,   1,   0,   1,   0,        57,      0,
     .                          0,    0,    0,    0,
     .            1,   0,  -2,   2,   0,       -58,      0,
     .                          0,    0,    0,    0,
     .            1,   1,   0,   0,   1,       -34,      0,
     .                         19,    0,    0,    0,
     .            2,   0,   0,   2,   0,        59,      0,
     .                          1,    0,    0,    0,
     .           -1,   0,   0,   2,   2,       -38,      0,
     .                         17,    0,    0,    0,
     .            0,   0,   0,  -1,   1,        33,      0,
     .                        -18,    0,    0,    0 /
      data IX14/   0,   1,   0,  -2,   1,       -33,      0,
     .                         18,    0,    0,    0,
     .           -1,   0,   2,  -2,   2,        36,      0,
     .                        -16,    0,    0,    0,
     .           -1,   1,   0,   0,   1,       -31,      0,
     .                         17,    0,    0,    0,
     .            1,   0,   2,   1,   2,        33,      0,
     .                        -14,    0,    0,    0,
     .            0,   0,   0,   4,   0,        48,      0,
     .                          1,    0,    0,    0,
     .            0,   0,   2,   1,   1,        27,      0,
     .                        -14,    0,    0,    0,
     .            1,   0,   0,  -2,   2,        32,      0,
     .                        -14,    0,    0,    0,
     .            1,   0,   2,  -1,   2,       -33,      0,
     .                         13,    0,    0,    0,
     .            1,  -1,   0,   2,   0,        48,      0,
     .                          0,    0,    0,    0,
     .           -1,   0,   2,   4,   1,       -26,      0,
     .                         13,    0,    0,    0 /
      data IX15/   0,   0,   2,   2,   0,        41,      0,
     .                          1,    0,    0,    0,
     .            1,   0,  -2,   0,   1,        27,      0,
     .                        -14,    0,    0,    0,
     .           -1,   0,   2,  -1,   1,       -23,      0,
     .                         14,    0,    0,    0,
     .            1,   1,   2,  -2,   1,        23,      0,
     .                        -12,    0,    0,    0,
     .            4,   0,   2,   0,   2,       -26,      0,
     .                         11,    0,    0,    0,
     .            0,   1,   2,   1,   2,       -24,      0,
     .                         10,    0,    0,    0,
     .            2,   0,   2,   0,   0,        36,      0,
     .                          1,    0,    0,    0,
     .            2,   1,   2,  -2,   2,        25,      0,
     .                        -10,    0,    0,    0,
     .            2,  -1,   0,   0,   0,        38,      0,
     .                          0,    0,    0,    0,
     .           -1,  -1,   0,   0,   1,        21,      0,
     .                        -12,    0,    0,    0 /
      data IX16/  -2,   0,   2,   2,   1,        22,      0,
     .                        -11,    0,    0,    0,
     .            0,   0,   0,   0,   3,       -22,      0,
     .                         10,    0,    0,    0,
     .            1,   0,   4,  -2,   2,        23,      0,
     .                         -9,    0,    0,    0,
     .            2,   0,   2,   2,   1,       -19,      0,
     .                         10,    0,    0,    0,
     .           -2,   0,   2,   4,   1,       -20,      0,
     .                         10,    0,    0,    0,
     .            0,   1,   0,   2,   1,        18,      0,
     .                         -9,    0,    0,    0,
     .            1,   0,   0,   1,   0,       -33,      0,
     .                          0,    0,    0,    0,
     .           -1,   0,   0,   4,   1,       -18,      0,
     .                          9,    0,    0,    0,
     .           -1,   0,   4,   0,   1,        19,      0,
     .                         -9,    0,    0,    0,
     .            0,   0,   2,  -3,   2,       -20,      0,
     .                          8,    0,    0,    0 /
      data IX17/   0,   0,   4,   0,   2,        19,      0,
     .                         -8,    0,    0,    0,
     .            2,   1,   0,   0,   0,       -28,      0,
     .                          0,    0,    0,    0,
     .            0,   0,   2,  -4,   1,       -16,      0,
     .                          9,    0,    0,    0,
     .           -1,  -1,   2,   4,   2,       -17,      0,
     .                          7,    0,    0,    0,
     .           -1,  -2,   0,   2,   0,        27,      0,
     .                          0,    0,    0,    0,
     .            0,   0,   0,   4,   1,       -16,      0,
     .                          7,    0,    0,    0,
     .            0,  -1,   0,   2,   1,       -14,      0,
     .                          7,    0,    0,    0,
     .            1,   0,   2,   4,   2,       -16,      0,
     .                          7,    0,    0,    0,
     .           -2,   0,   0,   2,   2,        18,      0,
     .                         -8,    0,    0,    0,
     .           -2,   2,   0,   2,   0,       -22,      0,
     .                          0,    0,    0,    0 /
      data IX18/  -2,  -1,   2,   0,   1,         9,      0,
     .                         -5,    0,    0,    0,
     .           -3,   0,   0,   0,   1,       -14,      0,
     .                          7,    0,    0,    0,
     .            0,   0,   2,   0,   3,        20,      0,
     .                          0,    0,    0,    0,
     .            0,   0,   2,   4,   1,       -12,      0,
     .                          6,    0,    0,    0,
     .            0,   0,   4,  -2,   1,        12,      0,
     .                         -7,    0,    0,    0,
     .            0,  -2,   0,   2,   0,        21,      0,
     .                          0,    0,    0,    0,
     .            1,   0,   0,  -1,   1,        17,      0,
     .                         -5,    0,   -3,    1,
     .            1,   1,   2,   2,   2,        15,      0,
     .                         -6,    0,    0,    0,
     .            3,   0,   2,  -2,   1,        12,      0,
     .                         -7,    0,    0,    0,
     .           -1,  -1,   2,   0,   2,       -16,      0,
     .                          6,    0,    0,    0 /
      data IX19/  -2,  -1,   0,   2,   1,       -13,      0,
     .                          7,    0,    0,    0,
     .            0,   0,   0,  -2,   2,        13,      0,
     .                         -5,    0,    0,    0,
     .            0,  -2,   2,   2,   2,       -13,      0,
     .                          5,    0,    0,    0,
     .            1,   0,   0,  -4,   1,       -12,      0,
     .                          6,    0,    0,    0,
     .           -1,   1,   0,   2,   1,       -10,      0,
     .                          6,    0,    0,    0,
     .           -2,   0,   0,   4,   1,        11,      0,
     .                         -6,    0,    0,    0,
     .            0,   0,   2,  -1,   1,       -10,      0,
     .                          5,    0,    0,    0,
     .            0,   2,   0,   0,   1,        -9,      0,
     .                          5,    0,    0,    0,
     .            0,   2,   2,  -2,   1,         8,      0,
     .                         -5,    0,    0,    0,
     .            2,   0,   0,   2,   1,        -9,      0,
     .                          5,    0,    0,    0 /
      data IX20/   2,   0,   0,  -4,   1,       -11,      0,
     .                          5,    0,    0,    0,
     .            2,   0,   2,  -4,   1,        10,      0,
     .                         -5,    0,    0,    0,
     .           -1,   0,  -2,   0,   1,       -10,      0,
     .                          5,    0,    0,    0,
     .           -1,   1,   2,   0,   1,         9,      0,
     .                         -5,    0,    0,    0,
     .           -1,   1,   2,  -2,   1,       -11,      0,
     .                          5,    0,    0,    0,
     .           -1,  -1,   0,   4,   0,        15,      0,
     .                          0,    0,    0,    0,
     .           -3,   0,   0,   4,   0,        16,      0,
     .                          0,    0,    0,    0,
     .            3,   0,   2,   2,   2,       -14,      0,
     .                          0,    0,    0,    0,
     .           -2,   1,   0,   2,   0,         9,      0,
     .                          1,    0,   -1,    0,
     .            0,   2,  -2,   2,   0,        -9,      0,
     .                          0,    0,    0,    0 /
      data IX21/   0,  -1,   2,   4,   2,        -9,      0,
     .                          0,    0,    0,    0,
     .            0,  -1,   2,  -1,   2,         9,      0,
     .                          0,    0,    0,    0,
     .            1,   1,   0,   2,   0,       -10,      0,
     .                          0,    0,    0,    0,
     .            2,   0,   0,  -2,   2,       -11,      0,
     .                          0,    0,    0,    0,
     .            2,  -1,   2,   2,   2,        -9,      0,
     .                          0,    0,    0,    0,
     .            4,   0,   0,   0,   0,         9,      0,
     .                          0,    0,    0,    0,
     .            4,   0,   2,  -2,   2,        12,      0,
     .                          0,    0,    0,    0,
     .           -1,   0,   0,   3,   0,       -10,      0,
     .                          0,    0,    0,    0,
     .           -1,   0,   4,  -2,   2,        -9,      0,
     .                          0,    0,    0,    0,
     .           -1,  -2,   2,   2,   2,        -9,      0,
     .                          0,    0,    0,    0 /
      data IX22/  -2,  -1,   0,   4,   0,        12,      0,
     .                          0,    0,    0,    0,
     .           -2,  -1,   2,   4,   2,       -12,      0,
     .                          0,    0,    0,    0,
     .            0,   1,   2,   2,   1,         7,      0,
     .                          0,    0,    0,    0,
     .            0,   2,   2,   0,   2,         7,      0,
     .                          0,    0,    0,    0,
     .            0,  -2,   2,   0,   2,        -8,      0,
     .                          0,    0,    0,    0,
     .            1,   0,   0,   4,   0,         8,      0,
     .                          0,    0,    0,    0,
     .            1,   0,   2,   2,   0,         8,      0,
     .                          0,    0,    0,    0,
     .            1,   0,   2,  -4,   2,         7,      0,
     .                          0,    0,    0,    0,
     .            1,  -1,   2,   2,   1,        -8,      0,
     .                          0,    0,    0,    0,
     .            1,  -1,   2,  -2,   2,        -7,      0,
     .                          0,    0,    0,    0 /
      data IX23/   1,  -2,   0,   0,   0,         8,      0,
     .                          0,    0,    0,    0,
     .            2,   0,   0,   0,   2,        -8,      0,
     .                          0,    0,    0,    0,
     .            2,   1,   0,  -2,   1,         8,      0,
     .                          0,    0,    0,    0,
     .            3,   0,   0,   0,   1,         7,      0,
     .                          0,    0,    0,    0,
     .           -1,   0,   2,   1,   2,         8,      0,
     .                          0,    0,    0,    0,
     .           -1,   0,   2,   3,   2,         8,      0,
     .                          0,    0,    0,    0,
     .           -1,   0,  -2,   4,   0,        -7,      0,
     .                          0,    0,    0,    0,
     .           -1,   1,   2,   2,   1,         7,      0,
     .                          0,    0,    0,    0,
     .           -1,   2,   0,   2,   0,        -8,      0,
     .                          0,    0,    0,    0,
     .           -1,  -1,   2,  -1,   1,         7,      0,
     .                          0,    0,    0,    0 /
      data IX24/  -2,   0,   2,  -2,   1,        -8,      0,
     .                          0,    0,    0,    0,
     .           -2,   0,   4,   0,   2,        -7,      0,
     .                          0,    0,    0,    0,
     .           -2,   0,  -2,   2,   0,         8,      0,
     .                          0,    0,    0,    0,
     .           -2,   1,   2,   0,   1,         9,      0,
     .                          0,    0,    0,    0,
     .           -3,   0,   2,   0,   1,        -8,      0,
     .                          0,    0,    0,    0,
     .            0,   1,   0,   1,   1,         5,      0,
     .                          0,    0,    0,    0,
     .            0,  -1,   0,   4,   0,         6,      0,
     .                          0,    0,    0,    0,
     .            0,  -1,   0,  -2,   1,         5,      0,
     .                          0,    0,    0,    0,
     .            0,  -2,   0,   0,   1,        -6,      0,
     .                          0,    0,    0,    0,
     .            1,   0,   2,   1,   1,         5,      0,
     .                          0,    0,    0,    0 /
      data IX25/   1,   0,   2,  -3,   2,        -6,      0,
     .                          0,    0,    0,    0,
     .            1,   0,  -2,   1,   0,        -7,      0,
     .                          0,    0,    0,    0,
     .            1,   1,   0,   1,   0,         5,      0,
     .                          0,    0,    0,    0,
     .            1,  -1,   0,  -2,   1,         6,      0,
     .                          0,    0,    0,    0,
     .            2,   0,   2,  -1,   2,        -6,      0,
     .                          0,    0,    0,    0,
     .            2,   1,   2,   0,   1,         5,      0,
     .                          0,    0,    0,    0,
     .            2,  -1,   2,   0,   1,        -6,      0,
     .                          0,    0,    0,    0,
     .            2,  -1,   2,  -2,   2,         5,      0,
     .                          0,    0,    0,    0,
     .            3,   0,   0,   2,   0,         5,      0,
     .                          0,    0,    0,    0,
     .            3,  -1,   2,   0,   2,        -5,      0,
     .                          0,    0,    0,    0 /
      data IX26/  -1,  -1,   2,   0,   1,        -6,      0,
     .                          0,    0,    0,    0,
     .           -2,   0,   0,   0,   2,         6,      0,
     .                          0,    0,    0,    0,
     .           -2,   0,   0,   3,   0,        -5,      0,
     .                          0,    0,    0,    0,
     .           -2,   0,   0,  -2,   1,        -5,      0,
     .                          0,    0,    0,    0,
     .           -2,   0,   2,   2,   0,        -6,      0,
     .                          0,    0,    0,    0,
     .           -2,  -1,   2,   0,   0,        -5,      0,
     .                          0,    0,    0,    0,
     .           -2,  -1,   2,   2,   2,         6,      0,
     .                          0,    0,    0,    0,
     .            0,   0,   1,   0,   0,         0,      0,
     .                          0,    0,    8,    0,
     .            0,   0,   1,   0,   1,         0,      0,
     .                          0,    0,  -16,  -14,
     .           -1,   0,   1,   0,   0,         0,      0,
     .                          0,    0,   33,    0 /
      data IX27/  -1,   0,   1,   0,   1,         0,      0,
     .                          0,    0, -105,  -89,
     .           -1,   0,   1,   0,   2,         0,      0,
     .                          0,    0,   36,   18,
     .           -1,   0,   1,   0,   3,         0,      0,
     .                          0,    0,   -6,    0  /
     
*           RFCN Mode   430.07 d 1979/ 1/ 1-1984/ 1/ 1
* Amplitudes:         .232     -.215      .000      .000mas
* Coefficients:      -.584     -.542     -.232      .215 mas
*           RFCN Mode   430.07 d 1984/ 1/ 1-1986/ 1/ 1
* Amplitudes:         .090     -.245      .000      .000mas
* Coefficients:      -.226     -.617     -.090      .245 mas
*           RFCN Mode   430.07 d 1986/ 1/ 1-1988/ 1/ 1
* Amplitudes:         .189     -.217      .000      .000mas
* Coefficients:      -.475     -.546     -.189      .217 mas
*           RFCN Mode   430.07 d 1988/ 1/ 1-1990/ 1/ 1
* Amplitudes:         .101     -.173      .000      .000mas
* Coefficients:      -.254     -.436     -.101      .173 mas
*           RFCN Mode   430.07 d 1990/ 1/ 1-1992/ 1/ 1
* Amplitudes:        -.025     -.170      .000      .000mas
* Coefficients:       .063     -.426      .025      .170 mas
*           RFCN Mode   430.07 d 1992/ 1/ 1-1994/ 1/ 1
* Amplitudes:        -.048     -.118      .000      .000mas
* Coefficients:       .120     -.296      .048      .118 mas
*           RFCN Mode   430.07 d 1994/ 1/ 1-1996/ 1/ 1
* Amplitudes:        -.004     -.076      .000      .000mas
* Coefficients:       .010     -.191      .004      .076 mas
     

****  Initialize the values and sum over the series
 
      dpsi_lsu = 0.0d0
      deps_lsu = 0.0d0
 
      cent = (epoch-DJ2000) / 36525.d0
 
      do 200 i = num_ls, 1, -1
 
*         Sum the mulitpliers by the arguments to the argument of
*         nutation
          arg = 0.d0
          do 150 j = 1,5
 
*            Sum into the argument for nutation.
             arg = arg + nutc_int(j,i)*ls_arg(j)
 150      continue
 
          arg = mod(arg, 2.d0*pi)
 
****      Now add contributions to dpsi and deps
          dpsi_lsu = dpsi_lsu +
     .               (nutc_int( 6,i)+ nutc_int(7,i)*cent)*sin(arg) +
     .                nutc_int(10,i)*cos(arg)
          deps_lsu = deps_lsu +
     .               (nutc_int( 8,i)+ nutc_int(9,i)*cent)*cos(arg) +
     .                nutc_int(11,i)*sin(arg)
 
 200  continue
 
*     Convert values from micro-arc-sec to mill-arc-second
      dpsi_ls = dpsi_lsu * 1.d-3
      deps_ls = deps_lsu * 1.d-3
 
****  Thats all
      return
      end
 
CTITLE OUT_PLAN_NUT
 
      subroutine out_plan_nut   

*     Routine to write the planetary contribution to the nutations
*     to stdout.
*     

* USAGE:
*     call out_plan_nut

* APPOXIMATIONS: The Oppolzer terms have not been added (should be
*                < 0.005 mas), and
*                Contributions from a non-rigid Earth have not been
*                computed.  For many of these terms the contribution
*                arises from the perturbation of the Earth's orbit and
*                therefore there will be not deformation effects.

 
* PASSED VARIABLES
*   NONE

* LOCAL VARIABLES
 
*   epoch       - Dummy Julian date (used in call to plan_angles).
*   plan_arg(10) - Values of the planetary arguments (Lve, Le,
*                 Lma,  LJ, Lsa, pa, D, F, lm, Om) (rads) 
*                 (same order as KS1990)
*   plan_rat(10) - Rates of changes of the planetary arguments
*                  (rad/year).  Used to get periods for the
*                 terms.
*   dpsi, deps   - Dummy returns for nutations in long and oblquity.
 
      real*8 epoch, plan_arg(10), plan_rat(10), dpsi, deps

****  Set the epoch to a dummy value

      epoch = 2 400 000.5d0 
 
***** Get the fundamental arguments at this epoch
 
      call plan_angles( epoch, plan_arg, plan_rat)

*     Now compute the contributions of the planetery ntations by 
*     summing over the series.

      call eval_plan_nut( plan_arg, plan_rat, dpsi, deps, 'YES' )
 
****  Thats all
      return
      end
 

CTITLE PLAN_NUT
 
      subroutine plan_nut( jd, dpsi, deps )   

*     Routine to compute the planetary contribution to the nutations.
*     Coefficents from Tables XIV to XIX of Kinoshita, H. and J. Souchay,
*     Nutations for the rigid Earth, Celes. Mech. and Dynam. Astron,
*
* NEW Version based on: Corrections and new developments in rigid Earth
*     nutation theory: Lunisolar influence including indirect planetary
*     effects, J. Souchay and H. Kinioshita, Astron. and Astrophys., 1995.
* (Version here based on data files: KSRE95_FIG_PSI.DAT KSRE95_FIG_EPS.DAT
*  and generated with ks_plan.f)
* MOD for KSV_1996_2: Corrected rate of change argument for lmc in 
*     plan_angles subroutine (see comments in routine) 
* MOD for KSV_1996_3: Replaced all planetary arguments with values from 
*     Simon et al., 1994.   



* USAGE:
*     call plan_nut( jd, dpsi, deps )   
*     where <jd>    is a full julian date with fractional part
*                   of the day added (REAL*8 INPUT)
*     and <dpsi> and <deps> are the contributions to the nutations
*                   in longitude and obliquity in milliarcsec.
*                   (REAL*8 OUTPUT)

* RESTRICTIONS: if <jd> is less than 2000000.0 this routine
*               assumes an MJD has been passed and the time
*               used will be converted to JD.  A warning 
*               message will be printed.
* APPOXIMATIONS: The Oppolzer terms have not been added (should be
*                < 0.005 mas), and
*                Contributions from a non-rigid Earth have not been
*                computed.  For many of these terms the contribution
*                arises from the perturbation of the Earth's orbit and
*                therefore there will be not deformation effects.

 
* PASSED VARIABLES
*
* INPUT Values
* jd     - Time at which value needed. (jd + fraction of day) 

* OUTPUT Values
* dpsi   - Contribution to the nutation in longitude (mas).  Should
*          be added to standard nutation in longitude values.
* deps   - Contribution to the nutation in obliquity (mas).  Should
*          be added to standard nutation in obliquity values.

      real*8 jd, dpsi, deps     

* LOCAL VARIABLES
 
*   epoch       - Julian date (jd passed in unless the JD 
*                 appears to be an MJD in which case it is 
*                 converted to JD (2 400 000.5d0 added) 
*   plan_arg(10) - Values of the planetary arguments (Lve, Le,
*                 Lma,  LJ, Lsa, pa, D, F, lm, Om) (rads) 
*                 (same order as KS1990)
*   plan_rat(10) - Rates of changes of the planetary arguments
*                  (rad/year).  Used to get periods for the
*                 terms.
 
      real*8 epoch, plan_arg(10), plan_rat(10)

***** Check to make sure user passed JD and not MJD.  Correct
*     problem and warn the user.
      if( jd .lt.2 000 000.0d0  ) then
          write(*,100) jd    
 100      format('**WARNING** MJD apparently passed to SD_COMP',
     .          ' Value (',F10.2,') converted to JD')
          epoch = jd + 2 400 000.5d0 
      else
          epoch = jd
      end if
 
***** Get the fundamental arguments at this epoch
 
      call plan_angles( epoch, plan_arg, plan_rat )

*     Now compute the contributions of the planetery ntations by 
*     summing over the series.

      call eval_plan_nut( plan_arg, plan_rat, dpsi, deps ,'NO')
 
****  Thats all
      return
      end

CTITLE EVAL_PLAN_NUT

      subroutine eval_plan_nut( plan_arg, plan_rat, dpsi, deps, out ) 

*     Routine to compute the planetary nutations by summing over the
*     KS1990 coefficients.  The coefficients and their arguments are
*     saved here in are integers in micro-arc-seconds.  

* NOTE: plan_angles must be called before routine.

* PARAMETERS:

* num_plan  - Number of contributions to the planetary nutations

      integer*4 num_plan 
      parameter ( num_plan      = 112 )

      real*8 pi
      parameter ( pi            = 3.1415926535897932D0 )

* PASSED PARAMETERS:

* INPUT:
* plan_arg(10)  - Ten planetary arguments including pa as given
*                 (KS1990.) (rad)
* plan_rat(10)  - Rates of change of the arguments (used to get
*           periods for the planetary nutations). (rad/yr)
* out           - Set to YES to output the planetary nutations, otherwize
*                 the actual angles will be computed

* OUTPUT:
* dpsi, deps   - Contributions to nutations in longuitude and 
*                obliquity (mas).

      real*8 plan_arg(10), plan_rat(10), dpsi, deps

      character*(*) out

* LOCAL VARIABLES:

* IX01-IX12(14,10) - Integer values of the planetary arguments and 
*                    values (micro-arc-seconds for values)
      integer*4 IX01(14,10), IX02(14,10), IX03(14,10),
     .          IX04(14,10), IX05(14,10), IX06(14,10),
     .          IX07(14,10), IX08(14,10), IX09(14,10),
     .          IX10(14,10), IX11(14,10), IX12(14, 2)

      integer*4 Plan_int(14,num_plan)

*  i and j   - Counters to loop over the coeffients and the argumemts

      integer*4 i,j

*  arg       - Final summed argumemt for the nutations contributions (rads)
*  dargdt    - Rate of change of the argument (rads/yr)
*  period    - Period of the nutation in days.
*  amp       - Total Amplitude of the planetary nutation.  (To be used for
*              sorting size)

      real*8 arg, dargdt, period, amp


      equivalence (Plan_int(1,  1),IX01)
      equivalence (Plan_int(1, 11),IX02)
      equivalence (Plan_int(1, 21),IX03)
      equivalence (Plan_int(1, 31),IX04)
      equivalence (Plan_int(1, 41),IX05)
      equivalence (Plan_int(1, 51),IX06)
      equivalence (Plan_int(1, 61),IX07)
      equivalence (Plan_int(1, 71),IX08)
      equivalence (Plan_int(1, 81),IX09)
      equivalence (Plan_int(1, 91),IX10)
      equivalence (Plan_int(1,101),IX11)
      equivalence (Plan_int(1,111),IX12)

      data IX01/   0,   2,   0,  -2,   0,   0,   2,   0,  -2,   1,
     .                                           28,    0,    0,  -15,
     .            18, -16,   0,   0,   0,   0,   0,   0,  -1,   0,
     .                                           23,   10,    0,    0,
     .            -8,  12,   0,   0,   0,   0,  -1,   1,   0,   1,
     .                                          120,   60,   32,  -64,
     .             0,   0,   2,   0,   0,   0,   1,  -1,   0,   0,
     .                                           27,   -8,    0,    0,
     .             0,  -1,   2,   0,   0,   0,   0,   0,   0,   1,
     .                                          -46,  -43,  -23,   25,
     .             3,  -4,   0,   0,   0,   0,   1,   0,  -1,   0,
     .                                            0,   13,    0,    0,
     .             5,  -6,   0,   0,   0,   0,   2,  -2,   0,   0,
     .                                            0,   23,    0,    0,
     .             6,  -8,   0,   0,   0,   0,   2,   0,  -2,   0,
     .                                            8,   -2,    0,    0,
     .             0,   8, -15,   0,   0,   0,   0,   0,   0,   0,
     .                                            5,   -2,    0,    0,
     .             0,  -2,   0,   3,   0,   0,  -2,   0,   2,   1,
     .                                            9,   -2,    0,   -5 /

      data IX02/   0,   2,   0,  -3,   0,   0,   2,   0,  -2,   0,
     .                                          -35,   -6,    0,    0,
     .             0,   1,   0,  -1,   0,   0,   1,   0,  -1,   0,
     .                                           -5,    0,    0,    0,
     .             0,   1,   0,   1,   0,   0,   1,  -1,   0,   0,
     .                                           -2,   -8,    0,    0,
     .             0,   0,   0,   1,   0,   0,   0,   0,   0,   1,
     .                                            2,   -7,    0,    0,
     .             0,  -1,   0,   0,  -1,   0,  -1,   1,   0,   1,
     .                                           17,    8,    4,   -9,
     .             0,   0,   0,   0,   1,   0,   0,   0,   0,   0,
     .                                            1,    6,    0,    0,
     .             0,   0,   0,   0,   1,   1,   0,   0,   0,   0,
     .                                            5,    0,    0,    0,
     .             0,  -1,   0,   0,   1,   0,  -1,   1,   0,   1,
     .                                           -7,   -1,    0,    0,
     .             8, -13,   0,   0,   0,   0,   0,   0,   0,   1,
     .                                            5,    7,    0,    0,
     .            18, -16,   0,   0,   0,   0,   0,   0,  -1,   1,
     .                                           -7,   -3,    0,    0 /

      data IX03/   0,   0,   0,  -2,   5,   0,   0,   0,   0,   1,
     .                                            8,    2,    0,    0,
     .             0,  -4,   8,  -3,   0,   0,   0,   0,   0,   1,
     .                                            8,   30,   16,   -4,
     .             0,   4,  -8,   3,   0,   0,   0,   0,   0,   1,
     .                                           -8,   29,   16,    4,
     .             0,   0,   0,   2,  -5,   0,   0,   0,   0,   1,
     .                                           -7,    2,    0,    0,
     .             0,   2,   0,  -2,   0,   0,   2,   0,  -2,   0,
     .                                          -44,    0,    0,    0,
     .           -18,  16,   0,   0,   0,   0,   0,   0,   1,   1,
     .                                            6,   -3,    0,    0,
     .            -8,  13,   0,   0,   0,   0,   0,   0,   0,   1,
     .                                           -4,    6,    0,    0,
     .             0,   0,  -2,   0,   0,   0,  -1,   1,   0,   1,
     .                                           27,    8,    4,  -15,
     .             0,   1,  -2,   0,   0,   0,   0,   0,   0,   0,
     .                                          -46,   44,    0,    0,
     .             0,  -2,   2,   0,   0,   0,  -1,   1,   0,   1,
     .                                            0,   -5,    0,    0 /

      data IX04/   0,   0,   0,   0,   2,   1,   0,   0,   0,   0,
     .                                            5,   10,    6,   -3,
     .             0,  -1,   0,   0,   2,   0,  -1,   1,   0,   1,
     .                                           -5,  -11,   -6,    3,
     .             0,   0,   0,   0,   2,   2,   0,   0,   0,   0,
     .                                          -12,    0,    0,    5,
     .            -5,   6,   0,   0,   0,   0,  -2,   2,   0,   1,
     .                                           -2,  -44,  -23,    0,
     .             0,   2,   0,  -3,   0,   0,   2,   0,  -2,   1,
     .                                           -5,    0,    0,    0,
     .             0,  -1,   0,  -1,   0,   0,  -1,   1,   0,   1,
     .                                           -5,   19,   10,    3,
     .             0,   0,   0,   1,   0,  -1,   0,   0,   0,   0,
     .                                            2,    6,    0,    0,
     .             0,   0,   0,   1,   0,   0,   0,   0,   0,   0,
     .                                           -8,   25,    0,    0,
     .             0,   0,   0,   1,   0,   1,   0,   0,   0,   0,
     .                                            0,    5,    0,    0,
     .             0,  -1,   0,   1,   0,   0,  -1,   1,   0,   1,
     .                                            0,   -5,    0,    0 /

      data IX05/   3,  -3,   0,   0,   0,   0,   2,   0,  -2,   0,
     .                                          -14,    0,    0,    0,
     .             0,  -2,   0,   2,   0,   0,  -2,   0,   2,   1,
     .                                            5,    0,    0,    0,
     .             3,  -5,   0,   0,   0,   0,   0,   0,   0,   0,
     .                                          -22,    7,    0,    0,
     .             3,  -5,   0,   0,   0,  -1,   0,   0,   0,   0,
     .                                           -1,   -7,    3,   -1,
     .             3,  -5,   0,   0,   0,  -2,   0,   0,   0,   0,
     .                                          211,    0,    0,   96,
     .             0,   2,  -4,   0,   0,   0,   0,   0,   0,   0,
     .                                           -8,   14,    0,    0,
     .             0,   2,  -4,   0,   0,  -2,   0,   0,   0,   0,
     .                                            5,    0,    0,    0,
     .             5,  -8,   0,   0,   0,  -2,   0,   0,   0,   0,
     .                                            0,  -26,   13,    0,
     .            -5,   7,   0,   0,   0,   0,  -1,   1,   0,   1,
     .                                           14,    3,    1,   -7,
     .             0,   0,   0,   2,   0,   1,   0,   0,   0,   0,
     .                                            4,   27,   12,   -2 /

      data IX06/   0,  -1,   0,   2,   0,   0,  -1,   1,   0,   1,
     .                                           -3,  -14,   -8,    1,
     .             0,   0,   0,   2,   0,   2,   0,   0,   0,   0,
     .                                         -116,    0,    0,   51,
     .            -3,   3,   0,   0,   0,   0,  -2,   2,   0,   1,
     .                                           12,    0,    0,   -6,
     .             0,  -2,   0,   2,   0,   0,  -2,   2,   0,   1,
     .                                            5,    0,    0,    0,
     .             2,  -3,   0,   0,   0,   0,   0,   0,   0,   0,
     .                                            0,   63,    0,    0,
     .             0,   0,   0,   3,   0,   2,   0,   0,   0,   0,
     .                                          -12,    0,    0,    5,
     .             1,  -2,   0,   0,   0,   0,   0,   0,   0,   0,
     .                                            0,   -9,    0,    0,
     .             0,   2,  -3,   0,   0,   0,   0,   0,   0,   0,
     .                                           -8,    5,    0,    0,
     .             0,   1,  -1,   0,   0,   0,   0,   0,   0,   0,
     .                                           -6,    0,    0,    0,
     .             4,  -7,   0,   0,   0,  -2,   0,   0,   0,   0,
     .                                            0,    6,    0,    0 /

      data IX07/   4,  -6,   0,   0,   0,  -2,   0,   0,   0,   0,
     .                                          -52,    0,    0,  -23,
     .             4,  -6,   0,   0,   0,  -1,   0,   0,   0,   0,
     .                                            0,    9,   -5,    1,
     .             1,  -1,   0,   0,   0,   0,   0,   0,   0,   0,
     .                                          153,    0,    0,    0,
     .             1,  -1,   0,   0,   0,   1,   0,   0,   0,   0,
     .                                            0,   -6,   -5,   -1,
     .             0,   1,   0,  -3,   0,  -2,   0,   0,   0,   0,
     .                                          -11,    0,    0,   -5,
     .             2,  -4,   0,   0,   0,  -1,   0,   0,   0,   0,
     .                                            0,   -8,    0,    0,
     .             2,  -4,   0,   0,   0,  -2,   0,   0,   0,   0,
     .                                           47,    0,    0,   20,
     .             0,   1,   0,  -2,   0,   0,   0,   0,   0,   0,
     .                                          -18,   27,    0,    0,
     .             0,   3,  -4,   0,   0,   0,   0,   0,   0,   0,
     .                                           -8,    5,    0,    0,
     .             3,  -4,   0,   0,   0,   0,   0,   0,   0,   0,
     .                                            0,   29,    0,    0 /

      data IX08/   0,   1,   0,  -1,   0,   0,   0,   0,   0,   0,
     .                                         -123,   -3,    0,    0,
     .             0,   2,  -2,   0,   0,   0,   0,   0,   0,   0,
     .                                          -38,    0,    0,    0,
     .             0,   1,   0,   0,  -1,   0,   0,   0,   0,   0,
     .                                           -8,    0,    0,    0,
     .             0,   0,   2,   0,   0,   2,   0,   0,   0,   0,
     .                                           -7,    0,    0,    0,
     .             0,   1,   0,   1,   0,   2,   0,   0,   0,   0,
     .                                          -25,    0,    0,   11,
     .             3,  -6,   0,   0,   0,  -2,   0,   0,   0,   0,
     .                                            0,   -5,    0,    0,
     .             5,  -7,   0,   0,   0,  -2,   0,   0,   0,   0,
     .                                          -21,    0,    0,   -9,
     .             0,   1,   0,   2,   0,   2,   0,   0,   0,   0,
     .                                           -3,   -5,    0,    0,
     .             2,  -2,   0,   0,   0,  -1,   0,   0,   0,   0,
     .                                            0,    5,    0,    0,
     .             2,  -2,   0,   0,   0,   0,   0,   0,   0,   0,
     .                                          -60,    0,    0,    0 /

      data IX09/   1,  -3,   0,   0,   0,  -2,   0,   0,   0,   0,
     .                                          -11,    0,    0,   -5,
     .             1,  -3,   0,   0,   0,  -1,   0,   0,   0,   0,
     .                                            0,   -5,    0,    0,
     .             0,   2,   0,  -3,   0,   0,   0,   0,   0,   0,
     .                                            8,    2,    0,    0,
     .             2,  -5,   0,   0,   0,  -2,   0,   0,   0,   0,
     .                                            0,  -13,    6,    0,
     .             6,  -8,   0,   0,   0,  -2,   0,   0,   0,   0,
     .                                          -12,    0,    0,   -5,
     .             0,   2,   0,  -2,   0,   0,   0,   0,   0,   0,
     .                                           39,    0,    0,    0,
     .             3,  -3,   0,   0,   0,   0,   0,   0,   0,   0,
     .                                           10,    0,    0,    0,
     .             3,  -3,   0,   0,   0,   2,   0,   0,   0,   0,
     .                                            5,   -2,    0,    0,
     .             0,   2,   0,  -1,   0,   2,   0,   0,   0,   0,
     .                                          -15,   -3,   -1,    7,
     .             0,   3,  -2,   0,   0,   2,   0,   0,   0,   0,
     .                                            8,   -7,    0,    0 /

      data IX10/   8, -15,   0,   0,   0,  -2,   0,   0,   0,   0,
     .                                           -6,  -10,    4,   -3,
     .             0,   6,  -8,   3,   0,   2,   0,   0,   0,   0,
     .                                           12,  -42,  -18,   -5,
     .             0,   2,   0,   0,   0,   2,   0,   0,   0,   0,
     .                                           -9,    0,    0,    0,
     .             0,   2,  -8,   3,   0,  -2,   0,   0,   0,   0,
     .                                           12,  -42,   18,    5,
     .             8, -11,   0,   0,   0,   2,   0,   0,   0,   0,
     .                                           -6,  -10,   -4,   13,
     .             0,   1,   2,   0,   0,   2,   0,   0,   0,   0,
     .                                           -8,   -7,    0,    0,
     .             0,   2,   0,   1,   0,   2,   0,   0,   0,   0,
     .                                           17,   -1,    0,   -7,
     .             3,  -7,   0,   0,   0,  -2,   0,   0,   0,   0,
     .                                            7,   -2,    0,    0,
     .             2,  -1,   0,   0,   0,   2,   0,   0,   0,   0,
     .                                            0,  -17,   -8,    0,
     .             7,  -9,   0,   0,   0,  -2,   0,   0,   0,   0,
     .                                           -7,    0,    0,    0 /

      data IX11/   4,  -4,   0,   0,   0,   0,   0,   0,   0,   0,
     .                                           11,    0,    0,    0,
     .             1,   1,   0,   0,   0,   2,   0,   0,   0,   0,
     .                                          -30,    0,    0,   13,
     .             0,   3,   0,  -2,   0,   2,   0,   0,   0,   0,
     .                                            7,   -9,   -4,   -3,
     .             3,  -2,   0,   0,   0,   2,   0,   0,   0,   0,
     .                                            0,  -11,   -5,    0,
     .             0,   3,   0,  -1,   0,   2,   0,   0,   0,   0,
     .                                           52,    2,    0,  -22,
     .             0,   4,  -2,   0,   0,   2,   0,   0,   0,   0,
     .                                           14,    0,    0,   -6,
     .             8, -10,   0,   0,   0,  -2,   0,   0,   0,   0,
     .                                           -5,    0,    0,    0,
     .             5,  -5,   0,   0,   0,   0,   0,   0,   0,   0,
     .                                            7,    0,    0,    0,
     .             2,   0,   0,   0,   0,   2,   0,   0,   0,   0,
     .                                           39,    0,    0,  -17,
     .             0,   4,   0,  -2,   0,   2,   0,   0,   0,   0,
     .                                          -18,    0,    0,    8 /

      data IX12/  18, -16,   0,   0,   0,   0,   0,   2,  -1,   2,
     .                                          -13,   -6,   -3,    5,
     .           -18,  16,   0,   0,   0,   0,   0,   2,   1,   2,
     .                                           13,   -6,   -3,   -5 /


*      Last line has only  2 values

****  Initialize the values and sum over the series

      dpsi = 0.0d0
      deps = 0.0d0

***** If we are to output then write out the header lines
      if( out.eq.'YES' ) then
          write(*,100)
 100      format('* Souchay and Kinoshita 1995 Planetary nutations',/,
     .           '* #  LVe  LEa  LMa LJu LSa  pa  D   F   lm  Om ',
     .           '  Period ',
     .           '       Longitude ',15x,'Obliquity',10x,' Total Amp',/,
     .   '*',40x,' (Solar  ',
     .           '      Sin        Cos',10x,' Sin        Cos ',/,
     .   '*',40x,'  days)  ',
     .           '     (mas)      (mas)',10x,'(mas)      (mas)',
     .        5x,'(mas)')
      end if

      do 200 i = num_plan, 1, -1
          
*         Sum the mulitpliers by the arguments to the argument of
*         nutation
          arg = 0.d0
          dargdt = 0.d0
          do 150 j = 1,10

*            Planetary values
             arg = arg + Plan_int(j,i)*plan_arg(j)
             dargdt = dargdt + Plan_int(j,i)*plan_rat(j)
 150      continue

          arg = mod(arg, 2.d0*pi)
          period = (2*pi/dargdt)*365.25d0

****      Now add contributions to dpsi and deps
          dpsi = dpsi + (Plan_int(11,i)*Sin(arg) +
     .                   Plan_int(12,i)*Cos(arg)) * 1.d-3
          deps = deps + (Plan_int(13,i)*Sin(arg) +
     .                   Plan_int(14,i)*Cos(arg)) * 1.d-3
          if( out.eq.'YES' ) then
              amp = sqrt( ((plan_int(11,i)**2+Plan_int(12,i)**2)*0.4**2+
     .                    plan_int(13,i)**2+Plan_int(14,i)**2)*1.d-6)
              write(*,175) i, (Plan_int(j,i), j=1,10), period,
     .                        (Plan_int(j,i)/1000.d0, j=11,14),
     .                         amp
 175          format(i4,10(1x,I3),1x,F8.1,1x,2(F7.3,1x,F7.3,2x),1x,F7.3)
          end if


 200  continue

****  Thats all
      return
      end

CTITLE 'PLAN_ANGLES'
 
      subroutine plan_angles( epoch, plan_arg, plan_rat )
 
 
*     Routine to compute of planetary arguments for planetary
*     nutation.  The longitudes of the major planets is computed
*     from Simon et al., 1994.
*     Units are kept the same as Simon et al., and converted during
*     caluclations.
* Reference:
*     Simon, J. L., Bretagnon, P., Chapront, J., Chapront-Touze, M.,
*         Francou, G., Laskar, J., 1994, "Numerical Expressions for 
*         Precession Formulae and Mean Elements for the Moon and
*         Planets," Astron. Astrophys., 282, pp. 663-683.


 
* PHYSICAL CONSTANTS NEEDED FOR SD_COMP

*   pi          - Define here to full precision
*   DJ2000      - Julian date of J2000
*   sec360      - number of seconds in 360 degreees.

      real*8 pi, DJ2000, sec360
 
      parameter ( pi            = 3.1415926535897932D0 )
      parameter ( DJ2000        = 2451545.d0           )
      parameter ( sec360        = 1296000.d0           )


* PASSED VARIABLES

* INPUT
* epoch  - Julian date plus fraction of a day
*
* OUTPUT
* plan_arg(10)  - Planetary arguments for longtitudes of
*           Venus, Earth, Mars, Jupiter, Saturn, pa, D, F, lm and
*           Om in order given (rads)
* plan_rat(10)  - Rates of change of the arguments (used to get
*           periods for the planetary nutations). (rad/yr)

      real*8 epoch, plan_arg(10), plan_rat(10)

* LOCAL VARIABLES 
*      cent         - Centuries since J2000.
* See comments in data statements about units for each type (from
* Simon et al., units) Final angles are returned in rads and rads/yr
*      vl           - Venus longitude (rads)
*      vlc(2)       - coefficients for computing vl
*      tl           - Earth longitude (rads)
*      tlc(2)       - coefficients for computing tl
*      ml           - Mars  longitude (rads)
*      mlc(2)       - coefficients for computing ml
*      jl           - Jupliter longitude (rads)
*      jlc(2)       - coefficients for computing jl
*      sl           - Saturn longitude (rads)
*      slc(2)       - coefficients for computing sl
*      pa           - pa (rads)
*      pac(2)       - Coefficients for computing pa from KS1990.
*                     (Values converted from rates and acceleration
*                      in 1000's year to centuries).
*      Dr           - Mean elongation of the Moon from the Sun (rads)
*      drc(2)       - Coefficients for computing dc
*      Fr           - Moon's mean longitude minus Om (rad)
*      Frc(2)       - Coefficients for computing Fr
*      lm           - Mean longitude of moon minus mean longitude
*                     of perigee (rads)
*      lmc(2)       - Coefficients for computing lm
*      Om           - Longitude of the ascending node of the moon's
*                     mean orbit on the elliptic (rad)
*      Omc          - Coefficients for computing Om
 
      real*8 cent, vl, vlc(2), tl, tlc(2), ml, mlc(2), jl, jlc(2), sl,
     .    slc(2), pa, pac(2), Dr, drc(2), Fr, frc(2), 
     .    lm, lmc(2), Om, Omc(2)
 
     
* MOD TAH 960625: Changed all the arguments to values from Simon et al.
*     1994.  The Mean elements referred to mean dynamical ecliptic and
*     equinox of data are from Table 5.9 of Simon et al. We use
*     the same units as in paper.

*                   (degrees)         ("/thousand years)
      data vlc  / 181.979 800 85d0, 2 106 691 666.319 89d0   /
      data tlc  / 100.466 456 83d0, 1 296 027 711.034 29d0   /
      data mlc  / 355.432 999 58d0,   689 101 069.330 69d0   /
      data jlc  /  34.351 518 74d0,   109 306 899.894 53d0   /
      data slc  /  50.077 444 30d0,    44 046 398.470 38d0   /

*     Pa is from Equation (6).  (IAU76 Masses)  
*                  "/thousand yrs  "/(thousand yrs)**2
      data pac  /    50 288.200d0, 111.202 2d0 /

*     Delaunay variables from Equation 3.5.b (values same as 
*     ls_angles).  NOTES: As in Simon et al., 1994 the angle
*     rates here are in Centuries not thousand years. 
*                        "                     "/Julian Cent.
      data drc  /  1 072 260.703 690d0,  1 602 961 601.2090d0 /
      data frc  /    335 779.526 232d0,  1 739 527 262.8478d0 /
      data lmc  /    485 868.249 036d0,  1 717 915 923.2178d0 /
      data omc  /    450 160.398 036d0,     -6 962 890.5431d0 /
 
***** Get number of Centuries since J2000
 
      cent = (epoch-DJ2000) / 36525.d0
 
*     Compute arguments 
*     For longitudes (degree-> seconds and time in thousand yrs)
*     Final result in radians.
      vl  = ((vlc(1)*3600.d0 + vlc(2)*cent/10.d0)/sec360)*2*pi
      tl  = ((tlc(1)*3600.d0 + tlc(2)*cent/10.d0)/sec360)*2*pi
      ml  = ((mlc(1)*3600.d0 + mlc(2)*cent/10.d0)/sec360)*2*pi
      jl  = ((jlc(1)*3600.d0 + jlc(2)*cent/10.d0)/sec360)*2*pi
      sl  = ((slc(1)*3600.d0 + slc(2)*cent/10.d0)/sec360)*2*pi
      
      pa  = ((pac(1)*cent/10 + pac(2)*(cent/10)**2)/sec360)*2*pi
      
      dr  = ((drc(1) + drc(2)*cent)/sec360)*2*pi
      fr  = ((frc(1) + frc(2)*cent)/sec360)*2*pi
      lm  = ((lmc(1) + lmc(2)*cent)/sec360)*2*pi
      om  = ((omc(1) + omc(2)*cent)/sec360)*2*pi

****  Now save the values
      plan_arg( 1) = vl
      plan_arg( 2) = tl 
      plan_arg( 3) = ml 
      plan_arg( 4) = jl 
      plan_arg( 5) = sl 
      plan_arg( 6) = pa 
      plan_arg( 7) = dr 
      plan_arg( 8) = fr 
      plan_arg( 9) = lm 
      plan_arg(10) = om 
      
***** Save the rates of change.
      plan_rat( 1) = ((vlc(2)/sec360)*2*pi)/1000.d0
      plan_rat( 2) = ((tlc(2)/sec360)*2*pi)/1000.d0
      plan_rat( 3) = ((mlc(2)/sec360)*2*pi)/1000.d0
      plan_rat( 4) = ((jlc(2)/sec360)*2*pi)/1000.d0
      plan_rat( 5) = ((slc(2)/sec360)*2*pi)/1000.d0
      
      plan_rat( 6) = (((pac(1) + 2*pac(2)*cent/10)/sec360)
     .                               *2*pi)/1000.d0
     
      plan_rat( 7) = (drc(2)/100.d0/sec360)*2*pi
      plan_rat( 8) = (frc(2)/100.d0/sec360)*2*pi
      plan_rat( 9) = (lmc(2)/100.d0/sec360)*2*pi
      plan_rat(10) = (omc(2)/100.d0/sec360)*2*pi

***** Thats all
      return
      end

 
CTITLE FCN_NUT
 
      subroutine fcn_nut ( jd, dpsi_fcn, deps_fcn )
 
*     Routine to compute the consttributions of the freely excited
*     FCN mode to the nutations in longitude and obliquity.
 
* USAGE:
*     call fcn_nut( jd, dpsi_fcn, deps_fcn )
*     where <jd>    is a full julian date with fractional part
*                   of the day added (REAL*8 INPUT)
*     and <dpsi_fcn> and <deps_fcn> are the contributions to the nutations
*                   in longitude and obliquity in milliarcsec.
*                   (REAL*8 OUTPUT)
 
* RESTRICTIONS: if <jd> is less than 2000000.0 this routine
*               assumes an MJD has been passed and the time
*               used will be converted to JD.  A warning
*               message will be printed.
 
* RESTRICTIONS: This term represents as free excitation mode and
*               therefore will change with time (in much the same
*               way that the Chandler Wobble changes).  The
*               frequency of the FCN used here is accurate, but
*               coefficients used will depend on time.  The values
*               are interpolated over the 1979-1995 interval with
*               the first or last values being used outside of these
*               times.
 
* PARAMETERS:
 
*   DJ2000      - Julian date of J2000
*   solar_to_sidereal   - Conversion from solar days to sidereal
*                 days.
*   num_fcn     - Number of FCN amplitudes (linear interpolation between
*                 values).
 
 
      real*8 DJ2000, pi, solar_to_sidereal

      integer*4 num_fcn
 
      parameter ( pi            = 3.1415926535897932D0 )
      parameter ( DJ2000        = 2451545.d0           )
      parameter ( solar_to_sidereal = 1.002737909d0 )
      parameter ( num_fcn       = 7 )
 
* PASSED VARIABLES
*
* INPUT Values
* jd     - Time at which value needed. (jd + fraction of day)
 
* OUTPUT Values
* dpsi_fcn   - Contribution to the nutation in longitude (mas).  Should
*          be added to standard nutation in longitude values.
* deps_fcn   - Contribution to the nutation in obliquity (mas).  Should
*          be added to standard nutation in obliquity values.
 
 
      real*8 jd, dpsi_fcn, deps_fcn
 
* LOCAL VARIABLES
 
*   epoch       - Julian date (jd passed in unless the JD
*                 appears to be an MJD in which case it is
*                 converted to JD (2 400 000.5d0 added)
 
*   fcn_freq    - Freqency of the FCN mode (cycles per sidreal day).
*   fcn_arg     - Argument for the fcn mode computed from J2000 (rad).
*   fcn_tabl(2, num_fcn)  - Amplitude for the fcn free exciation.  These
*                 are converted to nutation in longitude and obliquity.  (mas)
*   fcn_ampl(2) - Interpolated FCN amplitude (mas)
*   fcn_jd(num_fcn) - Starting epochs for the fcn amplitudes (JD) 
*   sine        - Sine of the mean obliquity of the ecliptic.  (A constant
*                 value can be used here since the changes are small for
*                 this constribtion i.e., between 1980 and 2000 the error
*                   in the nutation in longitude is only 0.05 micro-arc-sec
*   dt          - Time difference between epoch and tabular interval (days)
*   dt_tab      - Time difference in tables values.
*   dfcn_amp(2) - Change in FCN amplitude between tabular points (mas)
 
 
      real*8 epoch, fcn_freq, fcn_arg, fcn_ampl(2), 
     .       fcn_tabl(2,num_fcn), sine, fcn_jd(num_fcn),
     .       dt, dt_tab, dfcn_amp(2)

*   i           - A counter used in do loop to find the correct pair of
*                 amplitudes to use

      integer*4 i
 
      data  fcn_freq  /   -1.00231887299d0 /

*     Time dependent values 
      data  fcn_jd   / 2443874.50d0,  2445700.50d0, 2446431.50d0,
     .                 2447161.50d0,  2447892.50d0, 2448622.50d0,
     .                 2449353.50d0  / 

      data  fcn_tabl /  .232d0,  -.215d0,     .090d0,  -.245d0,
     .                  .189d0,  -.217d0,     .101d0,  -.173d0,
     .                 -.025d0,  -.170d0,    -.048d0,  -.118d0,
     .                 -.004d0,  -.076d0 / 

      data  sine      /   0.3977771203d0 /
 
***** Check to make sure user passed JD and not MJD.  Correct
*     problem and warn the user.
      if( jd .lt.2 000 000.0d0  ) then
          write(*,100) jd
 100      format('**WARNING** MJD apparently passed to SD_COMP',
     .          ' Value (',F10.2,') converted to JD')
          epoch = jd + 2 400 000.5d0
      else
          epoch = jd
      end if

****  Find out which table values we should use.
      if( epoch.le.fcn_jd(1) ) then
          fcn_ampl(1) = fcn_tabl(1,1)
          fcn_ampl(2) = fcn_tabl(2,1)
      else if( epoch.ge. fcn_jd(num_fcn) ) then
          fcn_ampl(1) = fcn_tabl(1,num_fcn)
          fcn_ampl(2) = fcn_tabl(2,num_fcn)
      else
          do 200 i = 1, num_fcn-1
             if( epoch.ge.fcn_jd(i) .and. epoch.lt.fcn_jd(i+1) ) then
                 dt = epoch - fcn_jd(i)
                 dt_tab = fcn_jd(i+1) - fcn_jd(i)
                 dfcn_amp(1) = fcn_tabl(1,i+1) - fcn_tabl(1,i)
                 dfcn_amp(2) = fcn_tabl(2,i+1) - fcn_tabl(2,i)
                 fcn_ampl(1) = fcn_tabl(1,i) + (dfcn_amp(1)/dt_tab)*dt
                 fcn_ampl(2) = fcn_tabl(2,i) + (dfcn_amp(2)/dt_tab)*dt
             end if
 200      continue
      end if

 
***** Get the argument for the FCN mode at this times
 
      fcn_arg = -2*pi*(1.d0+fcn_freq)*solar_to_sidereal*(epoch-DJ2000)
 
      dpsi_fcn = (-fcn_ampl(1)*sin(fcn_arg) + 
     .             fcn_ampl(2)*cos(fcn_arg))/sine

      deps_fcn = (-fcn_ampl(1)*cos(fcn_arg) - 
     .             fcn_ampl(2)*sin(fcn_arg))
 
***** Thats all
      return
      end
 
CTITLE PREC_NUT
 
      subroutine prec_nut( jd, dpsi_prec, deps_prec )
 
*     Routine to evaluate the corrections to the nutations in longitude
*     and obliquity due to the corrections to the IAU-1976 Luni-solar
*     precession constant and the secular rate of change of the obliquity
*     of the ecliptic.
 
* PARAMETERS:
 
*   DJ2000      - Julian date of J2000
 
 
      real*8 DJ2000
 
      parameter ( DJ2000        = 2451545.d0           )
 
* PASSED VARIABLES
*
* INPUT Values
* jd     - Time at which value needed. (jd + fraction of day)
 
* OUTPUT Values
* dpsi_prec   - Contribution to the nutation in longitude (mas).  Should
*          be added to standard nutation in longitude values. Value
*          valid only when the IAU-1976 precession constant used to
*          compute the transformation to mean system.
* deps_prec   - Contribution to the nutation in obliquity (mas).  Should
*          be added to standard nutation in obliquity values.
 
 
      real*8 jd, dpsi_prec, deps_prec
 
* LOCAL VARIABLES
 
*   epoch       - Julian date (jd passed in unless the JD
*                 appears to be an MJD in which case it is
*                 converted to JD (2 400 000.5d0 added)
*   cent        - Number of Julian centuries since J2000.0
*   DpsiDt      - Correction to precession constant as a
*                 linear rate of change of nutation in
*                 longitude. (arc-second/century)
*   DepsDt      - Correction to rate of change of oblquity
*                 (arc-second/century)
 
 
 
      real*8 epoch, cent, DpsiDt,  DepsDt
 
      data  DpsiDt  /  -0.2957d0  /
      data  DepsDt  /  -0.0227d0  /
 
***** Check to make sure user passed JD and not MJD.  Correct
*     problem and warn the user.
      if( jd .lt.2 000 000.0d0  ) then
          write(*,100) jd
 100      format('**WARNING** MJD apparently passed to SD_COMP',
     .          ' Value (',F10.2,') converted to JD')
          epoch = jd + 2 400 000.5d0
      else
          epoch = jd
      end if
 
****  Compute the number of centuries
 
      cent = (epoch - DJ2000)/36525.d0
 
      dpsi_prec = DpsiDt*cent*1000.d0
      deps_prec = DepsDt*cent*1000.d0
 
*     Thats all
      return
      end
 
 
CTITLE ls_iau76
 
      subroutine ls_iau76( epoch, ls_arg )
 
*     Routine to compute the value of the fundamental argument
*     for Brown's arguments.  Arguments based on the IERS
*     standards.

* PHYSICAL CONSTANTS
 
*   pi          - Define here to full precision
*   rad_to_deg  - Conversion from radians to degs.
*   DJ2000      - Julian date of J2000
*   sec360      - number of seconds in 360 degreees.
 
 
      real*8 pi, rad_to_deg, DJ2000, sec360
 
      parameter ( pi            = 3.1415926535897932D0 )
      parameter ( DJ2000        = 2451545.d0           )
      parameter ( sec360        = 1296000.d0           )
 
*     Computed quanities
      parameter ( rad_to_deg    = 180.d0   /pi         )
 
*-------------------------------------------------------------------
 
* PASSED VARIABLES
 
* INPUT
* epoch  - Julian date for arguments (jd + fraction of day, REAL*8)
 
* OUTPUT
* ls_arg(5) -  Brown's arguments (radians, REAL*8)
 
 
      real*8 epoch, ls_arg(5)
 
* LOCAL VARIABLES
*      cent             - Julian centuries to DJ2000.
*      el,eld           - Mean longitude of moon minus mean
*                       - longitude of moon's perigee (arcsec)
*      elc(5)           - Coefficients for computing el
*      elp,elpd         - Mean longitude of the sun minus mean
*                       - longitude of sun perigee (arcsec)
*      elpc(5)          - Coeffiecents for computing elp
*      f,fd             - Moon's mean longitude minus omega (sec)
*      fc(5)            - Coefficients for computing f
*      d,dd             - Mean elongation of the moon from the
*                       - sun (arcsec)
*      dc(5)            - coefficients for computing d
*      om,omd           - longitude of the ascending node of the
*                       - moon's mean orbit on the elliptic
*                       - measured from the mean equinox of date
*      omc(5)           - Coefficients for computing om.
 
 
      real*8 cent, el,eld, elc(5), elp, elpd, elpc(5),
     .    f,fd, fc(5), d,dd, dc(5), om,omd, omc(5)
 
****  DATA statements for the fundamental arguments.
 
      data elc    /     0.064d0,    31.310d0,    715922.633d0,
     .             485866.733d0,    1325.0d0 /
      data elpc   /    -0.012d0,    -0.577d0,   1292581.224d0,
     .            1287099.804d0,      99.0d0 /
      data fc     /     0.011d0,   -13.257d0,    295263.137d0,
     .             335778.877d0,    1342.0d0/
      data dc     /     0.019d0,    -6.891d0,    1105601.328d0,
     .            1072261.307d0,    1236.0d0/
      data omc    /     0.008d0,     7.455d0,    -482890.539d0,
     .             450160.280d0,      -5.0d0/
 
****  Get the number of centuries to current time
 
      cent = (epoch-dj2000) / 36525.d0
 
****  Compute angular arguments
      el = elc(1) * cent**3 + elc(2) * cent**2 + elc(3) * cent
     .          + elc(4) + mod( elc(5) * cent, 1.d0 ) * sec360
      el = mod( el, sec360 )
      eld = 3.d0 * elc(1) * cent**2 + 2.d0 * elc(2) * cent + elc(3)
     .      + elc(5) * sec360
c
      elp = elpc(1) * cent**3 + elpc(2) * cent**2 + elpc(3) * cent
     .     + elpc(4) + mod( elpc(5) * cent, 1.d0 ) * sec360
      elp = mod( elp, sec360 )
      elpd = 3.d0 * elpc(1) * cent**2 + 2.d0 * elpc(2) * cent + elpc(3)
     .       + elpc(5) * sec360
c
      f = fc(1) * cent**3 + fc(2) * cent**2 + fc(3) * cent
     .     + fc(4) + mod( fc(5) * cent, 1.d0 ) * sec360
      f = mod( f, sec360 )
      fd = 3.d0 * fc(1) * cent**2 + 2.d0 * fc(2) * cent + fc(3)
     .     + fc(5) * sec360
c
      d = dc(1) * cent**3 + dc(2) * cent**2 + dc(3) * cent
     .     + dc(4) + mod( dc(5) * cent, 1.d0 ) * sec360
      d = mod( d, sec360 )
      dd = 3.d0 * dc(1) * cent**2 + 2.d0 * dc(2) * cent + dc(3)
     .     + dc(5) * sec360
c
      om = omc(1) * cent**3 + omc(2) * cent**2 + omc(3) * cent
     .     + omc(4) + mod( omc(5) * cent, 1.d0 ) * sec360
      om = mod( om, sec360 )
      omd = 3.d0 * omc(1) * cent**2 + 2.d0 * omc(2) * cent + omc(3)
     .      + omc(5) * sec360
c
 
****  Now save the values.  Convert values from arcseconds to radians
 
      ls_arg(1) = el / (3600.d0*rad_to_deg)
      ls_arg(2) = elp/ (3600.d0*rad_to_deg)
      ls_arg(3) = f  / (3600.d0*rad_to_deg)
      ls_arg(4) = d  / (3600.d0*rad_to_deg)
      ls_arg(5) = om / (3600.d0*rad_to_deg)
 
***** Thats all
      return
      end
 
