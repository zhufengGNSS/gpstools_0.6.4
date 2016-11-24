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
% [version]: $Revision: 2 $ $Date: 06/07/08 1:01 $
%            Copyright(c) 2004-2006 by T.Takasu, all rights reserved
% [history]: 04/06/18   0.1  new
%-----------------------------------------------------------------------------*/
#include "mex.h"

/* unit vector of satellite fixed coordinate ---------------------------------*/
extern void SatFixed(const double *rsat, const double *rsun, double *ex,
					 double *ey, double *ez);

/* mex interface -------------------------------------------------------------*/
extern void mexFunction(int nargout, mxArray *argout[], int nargin,
					    const mxArray *argin[])
{
	double *rsat,*rsun,*ex,*ey,*ez;

	if (nargin<2||mxGetM(argin[0])<3||mxGetM(argin[1])<3)
		mexErrMsgTxt("argin error");
	if (nargout!=3) mexErrMsgTxt("argout error"); 
	rsat=mxGetPr(argin[0]);
	rsun=mxGetPr(argin[1]);
	argout[0]=mxCreateDoubleMatrix(3,1,mxREAL); ex=mxGetPr(argout[0]);
	argout[1]=mxCreateDoubleMatrix(3,1,mxREAL); ey=mxGetPr(argout[1]);
	argout[2]=mxCreateDoubleMatrix(3,1,mxREAL); ez=mxGetPr(argout[2]);
	SatFixed(rsat,rsun,ex,ey,ez);
}
