/*------------------------------------------------------------------------------
% [system] : GpsTools
% [module] : satellite antenna offset correction
% [func]   : calculate satellite antenna offset correction
% [argin]  : rsat = satellite postion(m) (eci)
%            rrcv = station position(m) (eci)
%            rsun = sun position(m) (eci)
%            apc  = satellite antenna offset(m) [x;y;z] (sallite fixed coordinate)
%           (apv) = satellite antenna pcv parameters 
%                 apv(n,1) = nadir angle (n-1) deg offset(m) (n=1:16)
% [argout] : apcs = satellite antenna offset(m)
% [note]   :
% [version]: $Revision: 8 $ $Date: 06/07/08 1:01 $
%            Copyright(c) 2004-2006 by T.Takasu, all rights reserved
% [history]: 04/06/01   0.1  new
%            04/06/26   0.2  add antenna pcv correction
%-----------------------------------------------------------------------------*/
#include "mex.h"

extern double SatApc(const double *rsat, const double *rrcv, const double *rsun,
					 const double *apc, const double *apv);

/* mex interface --------------------------------------------------------------*/
extern void mexFunction(int nargout, mxArray *argout[], int nargin,
					    const mxArray *argin[])
{
	double *rsat,*rrcv,*rsun,*apc,*apv=NULL,*apcs;
	
	if (nargin<4||mxGetM(argin[0])<3||mxGetM(argin[1])<3||mxGetM(argin[2])<3||
	    mxGetM(argin[3])<3)
		mexErrMsgTxt("argin error");
	
	rsat=mxGetPr(argin[0]);
	rrcv=mxGetPr(argin[1]);
	rsun=mxGetPr(argin[2]);
	apc =mxGetPr(argin[3]);
	if (nargin>4) {
		if (mxGetM(argin[3])<16) mexErrMsgTxt("argin error");
		apv=mxGetPr(argin[4]);
	}
	argout[0]=mxCreateDoubleMatrix(1,1,mxREAL); apcs=mxGetPr(argout[0]);
	
	*apcs=SatApc(rsat,rrcv,rsun,apc,apv);
}
