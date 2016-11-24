/*------------------------------------------------------------------------------
% [system] : GpsTools
% [module] : mjd to calender date/time
% [func]   : convert mjd to calender date/time
% [argin]  : mjd  = modified julan date (day)
%           (sec) = seconds of the day
% [argout] : dt   = date/time [year,month,day,hour,min,sec]
% [note]   : valid after 1582/10/10
% [version]: $Revision: 4 $ $Date: 06/07/08 1:01 $
%            Copyright(c) 2004-2006 by T.Takasu, all rights reserved
% [history]: 04/02/23  0.1  new
%-----------------------------------------------------------------------------*/
#include "mex.h"

extern void MjdToCal(const double *mjd, const double *sec, double *dt);

/* mex interface -------------------------------------------------------------*/
extern void mexFunction(int nargout, mxArray *argout[], int nargin,
					    const mxArray *argin[])
{
	double *mjd,*sec=NULL,*dt;

	if (nargin<1||mxGetM(argin[0])<1) mexErrMsgTxt("argin error");
	if (nargout>1) mexErrMsgTxt("argout error"); 
	mjd=mxGetPr(argin[0]);
	if (nargin>=2&&mxGetM(argin[1])>=1) sec=mxGetPr(argin[1]);
	argout[0]=mxCreateDoubleMatrix(1,6,mxREAL); dt=mxGetPr(argout[0]);
	MjdToCal(mjd,sec,dt);
}
