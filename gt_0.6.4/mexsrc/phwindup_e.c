/*------------------------------------------------------------------------------
% [system] : GpsTools
% [module] : phase windup effect correction
% [func]   : calculate phase windup effect correction
% [argin]  : rsat = satellite position(m) (eci)
%            rsun = sun position(m) (eci)
%            posr = station position (m) (ecef)
%            U    = eci to ecef transformation matrix
%            phwinp = previous correction(rad)
%           (rmode) = rcv attitude mode (0:fixed station,1:leo satellite)
% [argout] : phwin  = phase windup phase correction(rad)
% [note]   :
% [version]: $Revision: 16 $ $Date: 2008-12-12 15:49:30 +0900 (金, 12 12 2008) $
%            Copyright(c) 2004-2008 by T.Takasu, all rights reserved
% [history]: 04/06/18   0.1  new
%            04/08/12   0.2  add argin rmode
%            08/12/03   0.3  use 0 instead of nan as previous correction
%-----------------------------------------------------------------------------*/
#include "qtcmn.h"
#include "mex.h"

extern void EcefToGeod(const double *epos, double *gpos, double *E);

/* unit vector of satellite fixed coordinate in eci --------------------------*/
static void SatfToEcsf(const double *rsat, const double *rsun, double *ex,
					   double *ey)
{
	double esun[3],ez[3];
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
/* phase windup effect correction --------------------------------------------*/
extern double phwindup(const double *rsat, const double *rsun,
					   const double *posr, const double *U, const double *phwinp,
					   const double *rmode)
{
	double gpos[3],E[9],Ut[9],ex[3],ey[3],xr[3],yr[3],xs[3],ys[3];
	double k[3],kys[3],kyr[3],Ds[3],Dr[3],DsDr[3],phwin;
	int i;
	
	if ((int)rmode[0]==0) {
		EcefToGeod(posr,gpos,E);
		for (i=0;i<3;i++) {
			ex[i]=E[3*i+1];	/* north */
			ey[i]=-E[3*i];	/* west */
		}
	}
	else {
		CROSS(posr,posr+3,ey); NORMV(ey); /* cross-track */
		CROSS(ey,posr,ex); NORMV(ex); /* along-track */
	}
	Tr(U,Ut);
	Mv(Ut,ex,xr);
	Mv(Ut,ey,yr);
	SatfToEcsf(rsat,rsun,xs,ys);
	Mv(Ut,posr,k);
	for (i=0;i<3;i++) k[i]-=rsat[i];
	NORMV(k);
	CROSS(k,ys,kys);
	CROSS(k,yr,kyr);
	for (i=0;i<3;i++) {
		Ds[i]=xs[i]-k[i]*DOT(k,xs)-kys[i];
		Dr[i]=xr[i]-k[i]*DOT(k,xr)+kyr[i];
	}
	phwin=acos(DOT(Ds,Dr)/NORM(Ds)/NORM(Dr));
	CROSS(Ds,Dr,DsDr);
	if (DOT(k,DsDr)<0.0) phwin=-phwin;

	return 2.0*M_PI*floor((*phwinp-phwin)/2.0/M_PI+0.5)+phwin;
}
