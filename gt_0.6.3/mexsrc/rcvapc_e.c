/*------------------------------------------------------------------------------
% [system] : GpsTools
% [module] : receiver antenna offset
% [func]   : receiver antenna offset
% [argin]  : azel = azimuth/elevation angle(rad) [az,el]
%            apc1 = L1 phase center offset (m) [up;north;east]
%            apc2 = L2 phase center offset (m) [up;north;east]
%            ecc  = antenna deltas (m) [up;north;east]
%           (apv1) = L1 phase center variation (m) [el0;el5;el10;...;el90]
%           (apv2) = L2 phase center variation (m) [el0;el5;el10;...;el90]
% [argout] : apcr = receiver antenna offset (m) [apcr1,apcr2]
% [note]   :
% [version]: $Revision: 5 $ $Date: 06/07/08 1:01 $
%            Copyright(c) 2004-2006 by T.Takasu, all rights reserved
% [history]: 04/06/01   0.1  new
%            04/08/26   0.2  consider elevation range-out
%            04/09/01   0.3  inverse sign of antenna pcv value
%-----------------------------------------------------------------------------*/
#include "mex.h"
#include "qtcmn.h"

/* interpolate phase center variation ----------------------------------------*/
static double InterpApv(const double *apv, double el)
{
	int i=floor(el/5.0);
	double a=el/5.0-(double)i;
	return apv[i]*(1.0-a)+apv[i+1]*a;
}
/* receiver antenna offset ---------------------------------------------------*/
extern void RcvApc(const double *azel, const double *apc1, const double *apc2,
				   const double *ecc, const double *apv1, const double *apv2,
				   double *apcr)
{
	double el,es[3],apce1[3],apce2[3],cosel;
	el=MAX(azel[1],0.0);
	el=MIN(azel[1],M_PI/2.0);
	cosel=cos(el);
	es[0]=sin(azel[0])*cosel;	/* east */
	es[1]=cos(azel[0])*cosel;	/* north */
	es[2]=sin(el);				/* up */
	apce1[0]=apc1[2]+ecc[2];	/* east */
	apce1[1]=apc1[1]+ecc[1];	/* north */
	apce1[2]=apc1[0]+ecc[0];	/* up */
	apce2[0]=apc2[2]+ecc[2];	/* east */
	apce2[1]=apc2[1]+ecc[1];	/* north */
	apce2[2]=apc2[0]+ecc[0];	/* up */
	apcr[0]=DOT(apce1,es);
	apcr[1]=DOT(apce2,es);
	if (apv2!=NULL) {
		apcr[0]-=InterpApv(apv1,el*RAD2DEG);
		apcr[1]-=InterpApv(apv2,el*RAD2DEG);
	}
}
