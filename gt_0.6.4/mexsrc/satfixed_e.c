/*------------------------------------------------------------------------------
% [system] : GpsTools
% [module] : satellite fixed coordinate unit vectors
% [func]   : satellite fixed coordinate unit vectors
% [argin]  : rsat = satellite position (m)
%            rsun = sun position (m)
% [argout] : ex   = x unit vector (m)
%            ey   = y unit vector (m)
%            ez   = z unit vector (m)
% [note]   :
% [version]: $Revision: 12 $ $Date: 2008-11-25 10:02:15 +0900 (火, 25 11 2008) $
%            Copyright(c) 2004-2006 by T.Takasu, all rights reserved
% [history]: 04/06/18   0.1  new
%-----------------------------------------------------------------------------*/
#include "qtcmn.h"

/* unit vector of satellite fixed coordinate ---------------------------------*/
extern void SatFixed(const double *rsat, const double *rsun, double *ex,
					 double *ey, double *ez)
{
	double esun[3];
	esun[0]=rsun[0]-rsat[0];
	esun[1]=rsun[1]-rsat[1];
	esun[2]=rsun[2]-rsat[2];
	NORMV(esun);
	ez[0]=-rsat[0];
	ez[1]=-rsat[1];
	ez[2]=-rsat[2];
	NORMV(ez);
	CROSS(ez,esun,ey);
	NORMV(ey);
	CROSS(ey,ez,ex);
}
