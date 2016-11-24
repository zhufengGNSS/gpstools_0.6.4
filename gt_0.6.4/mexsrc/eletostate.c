/*------------------------------------------------------------------------------
% [system] : GpsTools
% [module] : satellite orbit element to position/velocity
% [func]   : transform satellite orbit elements to state(position/velocity)
%            calculate partial derivatives of state by orbit elements
% [argin]  : ele   = orbit elements(m,deg) [a,e,i,OMG,omg,M]
% [argout] : state = position/velocity(m,m/sec) [x;y;z;xdot;ydot:zdot]
%            (dsde) = partial derivatives of state by orbit elements(6x6)
% [note]   : error if a<=0,e<=0,1<=e
% [version]: $Revision: 12 $ $Date: 2008-11-25 10:02:15 +0900 (火, 25 11 2008) $
%            Copyright(c) 2004-2006 by T.Takasu, all rights reserved
% [history]: 04/01/05   0.1  new
*-----------------------------------------------------------------------------*/
#include "mex.h"

extern int EleToState(const double *ele, double *state, double *dsde);

/* mex interface -------------------------------------------------------------*/
extern void mexFunction(int nargout, mxArray *argout[], int nargin,
					    const mxArray *argin[])
{
	double *ele,*state,*dsde=NULL;

	if (nargin<1||mxGetN(argin[0])<6) mexErrMsgTxt("argin error");
	if (nargout>2) mexErrMsgTxt("argout error"); 
	ele=mxGetPr(argin[0]);
	argout[0]=mxCreateDoubleMatrix(6,1,mxREAL); state=mxGetPr(argout[0]);
	if (nargout>=2) {
		argout[1]=mxCreateDoubleMatrix(6,6,mxREAL); dsde=mxGetPr(argout[1]);
	}
	if (!EleToState(ele,state,dsde)) mexErrMsgTxt("orbit element error");
}
