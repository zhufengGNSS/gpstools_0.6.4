/*------------------------------------------------------------------------------
% [system] : GpsTools
% [module] : receiver antenna correction
% [func]   : receiver antenna correction
% [argin]  : azel = azimuth/elevation angle(rad) [az,el]
%            apc1 = L1 phase center offset (m) [up;north;east]
%            apc2 = L2 phase center offset (m) [up;north;east]
%            ecc  = antenna deltas (m) [up;north;east]
%           (apv1)= L1 phase center variation (m) (73x19)
%           (apv2)= L2 phase center variation (m) (73x19)
%                   apv?[i+j*73]=az(i),el(j) pcv (az=0:5:360deg,el=0:5:90deg)
% [argout] : apcr = receiver antenna correction (m) [apcr1,apcr2]
%                   (observed range = geometric range + apcr)
% [note]   :
% [version]: $Revision: 16 $ $Date: 2008-12-12 15:49:30 +0900 (金, 12 12 2008) $
%            Copyright(c) 2004-2008 by T.Takasu, all rights reserved
% [history]: 04/06/01   0.1  new
%            04/08/26   0.2  consider elevation range-out
%            04/09/01   0.3  invert sign of antenna pcv value
%            08/11/25   0.4  invert sign of argout apcr (gt_0.6.4)
%                            support pcv azimuth angle dependancy 
%-----------------------------------------------------------------------------*/
#include "mex.h"
#include "qtcmn.h"

#define NV		73				/* number of rows in apv matrix */

/* interpolate phase center variation ----------------------------------------*/
static double InterpApv(const double *apv, const double *azel)
{
	double a,e,az=azel[0]*RAD2DEG,el=azel[1]*RAD2DEG;
	int i,j;
	az-=floor(az/360.0)*360; /* 0<=az<360 */
	if (el<0.0) el=0.0; else if (el>89.999) el=89.999; /* 0<=el<90 */
	i=(int)floor(az/5.0);
	j=(int)floor(el/5.0);
	a=az/5.0-(double)i;
	e=el/5.0-(double)j;
	if (0<=i&&i<72&&0<=j&&j<18) {
		return apv[i+j*NV]*(1.0-a)*(1.0-e)+apv[i+1+j*NV]*a*(1.0-e)+
			   apv[i+(j+1)*NV]*(1.0-a)*e+apv[i+1+(j+1)*NV]*a*e;
	}
	return 0.0;
}
/* receiver antenna offset ---------------------------------------------------*/
extern void RcvApc(const double *azel, const double *apc1, const double *apc2,
				   const double *ecc, const double *apv1, const double *apv2,
				   double *apcr)
{
	double es[3],apce1[3],apce2[3],cosel;
	cosel=cos(azel[1]);
	es[0]=sin(azel[0])*cosel;	/* east */
	es[1]=cos(azel[0])*cosel;	/* north */
	es[2]=sin(azel[1]);			/* up */
	apce1[0]=apc1[2]+ecc[2];	/* east */
	apce1[1]=apc1[1]+ecc[1];	/* north */
	apce1[2]=apc1[0]+ecc[0];	/* up */
	apce2[0]=apc2[2]+ecc[2];	/* east */
	apce2[1]=apc2[1]+ecc[1];	/* north */
	apce2[2]=apc2[0]+ecc[0];	/* up */
	apcr[0]=-DOT(apce1,es);
	apcr[1]=-DOT(apce2,es);
	if (apv1!=NULL) {
		apcr[0]+=InterpApv(apv1,azel);
	}
	if (apv2!=NULL) {
		apcr[1]+=InterpApv(apv2,azel);
	}
}
