/*------------------------------------------------------------------------------
% [system] : GpsTools
% [module] : multipath bias correction
% [func]   : multipath bias correction
% [argin]  : azel = azimuth/elevation angle(rad) [az,el]
%            mpcc = multipath bias cooefficient Cnm
%            mpcs = multipath bias cooefficient Snm
% [argout] : mpr  = multipath bias (m)
% [note]   :
% [version]: $Revision: 5 $ $Date: 06/07/08 1:01 $
%            Copyright(c) 2004-2006 by T.Takasu, all rights reserved
% [history]: 04/10/19   0.1  new
%-----------------------------------------------------------------------------*/
#include "mex.h"

extern double RcvMpc(const double *azel, const double *mpcc, const double *mpcs,
					 int nmax);

/* mex interface -------------------------------------------------------------*/
extern void mexFunction(int nargout, mxArray *argout[], int nargin,
					    const mxArray *argin[])
{
	double *azel,*mpcc,*mpcs,*mpr;
	int nmax;
	
	if (nargin<3||mxGetN(argin[0])<2||mxGetM(argin[1])!=mxGetN(argin[1])||
		mxGetM(argin[2])!=mxGetN(argin[2])||mxGetM(argin[1])!=mxGetM(argin[2]))
		mexErrMsgTxt("argin error");
	
	azel=mxGetPr(argin[0]);
	mpcc=mxGetPr(argin[1]);
	mpcs=mxGetPr(argin[2]);
	nmax=mxGetM(argin[1])-1;
	argout[0]=mxCreateDoubleMatrix(1,1,mxREAL); mpr=mxGetPr(argout[0]);
	
	*mpr=RcvMpc(azel,mpcc,mpcs,nmax);
}
