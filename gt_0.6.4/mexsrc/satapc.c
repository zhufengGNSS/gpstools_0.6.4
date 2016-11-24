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
%            08/11/25   0.3  invert sign of apcs (gt_0.6.4)
%                            support L1/L2 antenna correction
%-----------------------------------------------------------------------------*/
#include "mex.h"

extern void SatApc(const double *rsat, const double *rrcv, const double *rsun,
				   const double *apc1, const double *apc2, const double *apv1,
				   const double *apv2, double *apcs);

/* mex interface --------------------------------------------------------------*/
extern void mexFunction(int nargout, mxArray *argout[], int nargin,
					    const mxArray *argin[])
{
	double *rsat,*rrcv,*rsun,*apc1,*apc2,*apv1=NULL,*apv2=NULL,*apcs;
	
	if (nargin<5||mxGetM(argin[0])<3||mxGetM(argin[1])<3||mxGetM(argin[2])<3||
	    mxGetM(argin[3])<3||mxGetM(argin[4])<3)
		mexErrMsgTxt("argin error");
	
	rsat=mxGetPr(argin[0]);
	rrcv=mxGetPr(argin[1]);
	rsun=mxGetPr(argin[2]);
	apc1=mxGetPr(argin[3]);
	apc2=mxGetPr(argin[4]);
	if (nargin>5) {
		if (mxGetM(argin[5])<16) mexErrMsgTxt("argin error");
		apv1=mxGetPr(argin[5]);
	}
	if (nargin>6) {
		if (mxGetM(argin[6])<16) mexErrMsgTxt("argin error");
		apv2=mxGetPr(argin[6]);
	}
	argout[0]=mxCreateDoubleMatrix(1,2,mxREAL); apcs=mxGetPr(argout[0]);
	
	SatApc(rsat,rrcv,rsun,apc1,apc2,apv1,apv2,apcs);
}
