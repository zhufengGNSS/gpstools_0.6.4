/*------------------------------------------------------------------------------
% [system] : GpsTools
% [module] : calender date/time to mjd
% [func]   : convert calender date/time to mjd
% [argin]  : dt   = date/time [year,month,day(,hour,min,sec)]
% [argout] : mjd  = modified julan date (day)
%           (sec) = seconds of the day
% [note]   : valid after 1582/10/10
% [version]: $Revision: 16 $ $Date: 2008-12-12 15:49:30 +0900 (金, 12 12 2008) $
%            Copyright(c) 2004-2008 by T.Takasu, all rights reserved
% [history]: 04/02/23  0.1  new
%            08/12/04  0.2  suppress warning
%-----------------------------------------------------------------------------*/
#include "mex.h"

extern void CalToMjd(const double *dt, double *mjd, double *sec, int len_dt);

/* mex interface -------------------------------------------------------------*/
extern void mexFunction(int nargout, mxArray *argout[], int nargin,
					    const mxArray *argin[])
{
	double *dt,*mjd,*sec=NULL;

	if (nargin<1||mxGetN(argin[0])<3) mexErrMsgTxt("argin error");
	if (nargout>2) mexErrMsgTxt("argout error"); 
	dt=mxGetPr(argin[0]);
	argout[0]=mxCreateDoubleMatrix(1,1,mxREAL); mjd=mxGetPr(argout[0]);
	if (nargout>1) {
		argout[1]=mxCreateDoubleMatrix(1,1,mxREAL); sec=mxGetPr(argout[1]);
	}
	CalToMjd(dt,mjd,sec,(int)mxGetN(argin[0]));
}
