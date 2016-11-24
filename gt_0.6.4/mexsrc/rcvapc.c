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

extern void RcvApc(const double *azel, const double *apc1, const double *apc2,
				   const double *ecc, const double *apv1, const double *apv2,
				   double *apcr);

/* mex interface -------------------------------------------------------------*/
extern void mexFunction(int nargout, mxArray *argout[], int nargin,
					    const mxArray *argin[])
{
	double *azel,*apc1,*apc2,*ecc,*apv1=NULL,*apv2=NULL,*apcr;
	
	if (nargin<4||mxGetN(argin[0])<2||mxGetM(argin[1])<3||mxGetM(argin[2])<3||
	    mxGetM(argin[3])<3)
		mexErrMsgTxt("argin error");
	
	azel=mxGetPr(argin[0]);
	apc1=mxGetPr(argin[1]);
	apc2=mxGetPr(argin[2]);
	ecc =mxGetPr(argin[3]);
	if (nargin>4&&mxGetM(argin[4])>=19) apv1=mxGetPr(argin[4]);
	if (nargin>5&&mxGetM(argin[5])>=19) apv2=mxGetPr(argin[5]);
	argout[0]=mxCreateDoubleMatrix(1,2,mxREAL); apcr=mxGetPr(argout[0]);
	
	RcvApc(azel,apc1,apc2,ecc,apv1,apv2,apcr);
}
