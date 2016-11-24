/*------------------------------------------------------------------------------
% [system] : GpsTools
% [module] : satellite state to orbit element
% [func]   : convert satellite position/velocity to orbit element
% [argin]  : state = position/velocity(m,m/sec) [x;y;z;xdot;ydot:zdot]
% [argout] : ele   = orbit element(m,deg) [a,e,i,OMG,omg,M]
% [note]   : tangental orbit element
% [version]: $Revision: 12 $ $Date: 2008-11-25 10:02:15 +0900 (火, 25 11 2008) $
%            Copyright(c) 2004-2006 by T.Takasu, all rights reserved
% [history]: 04/01/05   0.1  new
%----------------------------------------------------------------------------*/
#include "mex.h"
#include "qtcmn.h"

extern void StateToEle(const double *state, double *ele);

/* mex interface ------------------------------------------------------------*/
extern void mexFunction(int nargout, mxArray *argout[], int nargin,
					    const mxArray *argin[])
{
	double *state,*ele;

	if (nargin<1||mxGetM(argin[0])<6) mexErrMsgTxt("argin error");
	if (nargout>1) mexErrMsgTxt("argout error"); 
	state=mxGetPr(argin[0]);
	argout[0]=mxCreateDoubleMatrix(1,6,mxREAL); ele=mxGetPr(argout[0]);
	
	StateToEle(state,ele);
}
