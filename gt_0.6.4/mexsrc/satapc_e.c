/*------------------------------------------------------------------------------
% [system] : GpsTools
% [module] : satellite antenna correction
% [func]   : satellite antenna correction
% [argin]  : rsat = satellite postion(m) (eci)
%            rrcv = station position(m) (eci)
%            rsun = sun position(m) (eci)
%            apc1 = L1 phase center offset (m) [x;y;z]
%            apc2 = L2 phase center offset (m) [x;y;z]
%           (apv1)= L1 phase center variation (m) (1x16)
%           (apv2)= L2 phase center variation (m) (1x16)
%                   apv?[i]=nadir(i) pcv (nadir=0:1:15deg)
% [argout] : apcs = satellite antenna correction (m) [apcs1,apcs2]
%                   (observed range = geometric range + apcs)
% [note]   :
% [version]: $Revision: 16 $ $Date: 2008-12-12 15:49:30 +0900 (金, 12 12 2008) $
%            Copyright(c) 2004-2008 by T.Takasu, all rights reserved
% [history]: 04/06/01   0.1  new
%            04/06/26   0.2  add antenna pcv correction
%            08/11/25   0.3  invert sign of apcs (gt_0.6.4)
%                            support L1/L2 antenna correction
%                            fix bug on error-stop if isnan(rrcv)
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
	nadir=180.0-acos(DOT(rs,rr)/NORM(rs))*RAD2DEG; /* nadir angle (deg) */
	n=(int)floor(nadir);
	a=nadir-(double)n;
	return n>=15?apv[15]:(n>=0?apv[n]*(1.0-a)+apv[n+1]*a:0.0);
}
/* satellite antenna offset correction ---------------------------------------*/
extern void SatApc(const double *rsat, const double *rrcv, const double *rsun,
				   const double *apc1, const double *apc2, const double *apv1,
				   const double *apv2, double *apcs)
{
	double E[9],off1[3],off2[3],rr[3];
	
	SatfToEcsf(rsat,rsun,E);
	Mv(E,apc1,off1);
	Mv(E,apc2,off2);
	rr[0]=rrcv[0]-rsat[0];
	rr[1]=rrcv[1]-rsat[1];
	rr[2]=rrcv[2]-rsat[2];
	NORMV(rr);
	apcs[0]=-DOT(off1,rr);
	apcs[1]=-DOT(off2,rr);
	if (apv1!=NULL) {
		apcs[0]+=SatApv(apv1,rsat,rr);
	}
	if (apv2!=NULL) {
		apcs[1]+=SatApv(apv2,rsat,rr);
	}
}
