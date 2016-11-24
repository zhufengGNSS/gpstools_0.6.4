%-------------------------------------------------------------------------------
% [system] : GpsTools
% [module] : read grib gpv data file
% [func]   : read grib gpv data file
% [argin]  : file   = grib file path
% [argout] : data   = gpv data (struct array)
%                data(n).pprm = product parameter struct
%                  cid   : center id
%                  sid   : subcenter id
%                  pid   : process id
%                  gid   : grid id
%                          21-64 : international exchange grids (see GRIB specs)
%                          255 : non-defined grids
%                  param : parameter and units
%                            2 : mean sea level pressure(Pa)
%                            7 : geopotentiol height(m)
%                           11 : temperture(K)
%                           33 : wind u component(m/s)
%                           34 : wind v component(m/s)
%                           39 : vertical velocity(Pa/s)
%                           52 : relative humidity(%)
%                           81 : land ratio(1:land,0:sea,rate)
%                  level : layer level
%                          [  1,0] : surface
%                          [100,p] : pressure in p hPa
%                          [102,0] : mean sea level
%                  time  : initial time of forecast (UTC)
%                          [year,month,day,hour,min,sec]
%                  tunit : forecast time unit
%                            0 : min
%                            1 : hour
%                            2 : day
%                  p1    : period of time in tunit (0:analysis/initialized)
%                  p2    : period of time or time interval in tunit
%                data(n).gprm = grid parameter struct
%                  nv    : number of vertical coordinate parameters
%                  pv    : location of vertical coordinate parameters
%                  type  : data representation type
%                            0 : latitude/longitude grid (eq-cylindrical)
%                            1 : mercator projection
%                            3 : lambert conformal projection
%                            5 : polar stereographic projection
%                           13 : oblique lambert conformal projection
%                           50 : spherical harmonic coefficients
%                           90 : space view of orthgraphic grid
%                  nx,ny : no. of points along lat/lon or x/y-axis
%                  rcflg : resolution and component flag
%                           rcflg(1) : direction increment given
%                           rcflg(2) : earth radius
%                               0 : sphere re = 6367.47km
%                               1 : oblate spheriod re=6378.16km f=1/297.0
%                           rcflg(3) : uv component of vector
%                               0 : relative to east/north
%                               1 : relative to defined grid x/y
%                  smode : scanning mode
%                           0 : points scan in + direction
%                           1 : points scan in - direction
%                  lat1,lon1 : lat/lon of first grid (deg or m)
%                  lat2,lon2 : lat/lon of last grid (deg) (type=0)
%                  dx,dy : lat/lon or x/y increment (deg or m)
%                  lov   : orientation of grid (y-parallel lon.) (deg) (type!=0)
%                  lati1,lati2 : lats which secant cone cuts the sphere (type!=0)
%                  latp,lonp : lat/lon of southern pole (type!=0)
%                data(n).data = grid data
% [note]   :
% [version]: $Revision: 12 $ $Date: 2008-11-25 10:02:15 +0900 (火, 25 11 2008) $
%            Copyright(c) 2004-2006 by T.Takasu, all rights reserved
% [history]: 04/03/08   0.1  new
%-----------------------------------------------------------------------------*/

% (mex function)

