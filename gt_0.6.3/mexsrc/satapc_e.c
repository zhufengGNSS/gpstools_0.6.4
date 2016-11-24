/*------------------------------------------------------------------------------
% [system] : GpsTools
% [module] : satellite antenna offset
% [func]   : calculate satellite antenna offset correction
% [argin]  : rsat = satellite postion(m) (eci)
%            rrcv = station position(m) (eci)
%            rsun = sun position(m) (eci)
%            apc  = satellite antenna offset(m) [x;y;z] (sallite fixed coordinate)
%           (apv) = satellite antenna pcv parameters
%                 apv(n,1) = nadir angle (n-1) deg offset(m) (n=1:16)
% [argout] : apcs = satellite antenna offset(m)
% [note]   :
% [version]: $Revision: 6 $ $Date: 06/07/08 1:01 $
%            Copyright(c) 2004-2006 by T.Takasu, all rights reserved
% [history]: 04/06/01   0.1  new
%            04/06/26   0.2  add antenna pcv correction
%-----------------------------------------------------------------------------*/
#include "qtcmn.h"

/* satellite fixed coordinate to eci transformation matrix -------------------*/
static void SatfToEcsf(const double *rsat, const double *rsun, double *E)
{
	double esun[3],*ex=E,*ey=E+3,*ez=E+6;
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
/* satellite antenna pcv correction ------------------------------------------*/
extern double SatApv(const double *apv, const double *rs, const double *rr)
{
	double nadir,a;
	int n;
	if (apv==NULL) return 0.0;
	nadir=180.0-acos(DOT(rs,rr)/NORM(rs))*RAD2DEG; /* nadir angle (deg) */
	n=(int)nadir;
	a=nadir-(double)n;
	return n>=15?apv[15]:apv[n]*(1.0-a)+apv[n+1]*a;
}
/* satellite antenna offset correction ---------------------------------------*/
extern double SatApc(const double *rsat, const double *rrcv, const double *rsun,
					 const double *apc, const double *apv)
{
	double E[9],off[3],rr[3];
	SatfToEcsf(rsat,rsun,E);
	Mv(E,apc,off);
	rr[0]=rrcv[0]-rsat[0];
	rr[1]=rrcv[1]-rsat[1];
	rr[2]=rrcv[2]-rsat[2];
	NORMV(rr);
	return DOT(off,rr)+SatApv(apv,rsat,rr);
}
