%-------------------------------------------------------------------------------
% [system] : GpsTools
% [module] : read VMF1 grid file
% [func]   : read VMF1 grid file
% [argin]  : file = filepath
% [argout] : data  = data(i,j) lats(i) lons(j) grid data (nxm)
%            lats  = latitudes  (deg) (1xn)
%            lons  = longitudes (deg) (1xm)
% [note]   : reference :
%            J.Boehm et al., Troposphere mapping functions for GPS and very long
%            baseline interferometry from European Centre for Medium-Range
%            Weather Forecasts operational analysis data, J. Geoph. Res.,
%            Vol. 111, B02406, 2006
%            vmf1 grid file can be downloaded from :
%            http://www.hg.tuwien.ac.at/~ecmwf1/
% [version]: $Revision: 3 $ $Date: 06/07/08 1:01 $
%            Copyright(c) 2004-2008 by T.Takasu, all rights reserved
% [history]: 08/12/06  0.1  new
%-------------------------------------------------------------------------------

% (mex function)

