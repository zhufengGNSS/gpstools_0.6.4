/*------------------------------------------------------------------------------
% [system] : GpsTools
% [module] : eci to satellite-fixed coordinate transformation matrix
% [func]   : calculate eci to satellite-fixed coordinate transformation matrix
% [argin]  : state = satellite position/velocity[x;y;z;vx;vy;vz](m) (eci)
% [argout] : E = eci to satellite-fixed coordinate transformation matrix(3x3)
%                [r_radial;r_along-track;r_cross-track]=E*r
% [note]   :
% [version]: $Revision: 12 $ $Date: 2008-11-25 10:02:15 +0900 (火, 25 11 2008) $
%            Copyright(c) 2004-2006 by T.Takasu, all rights reserved
% [history]: 04/02/08   0.1  new
%            04/08/11   0.2  separated from ecsftosatf_e.c
%-----------------------------------------------------------------------------*/
#include "mex.h"

extern void EcsfToSatf(const double *state, double *E);

/* mex interface -------------------------------------------------------------*/
extern void mexFunction(int nargout, mxArray *argout[], int nargin,
					    const mxArray *argin[])
{
	double *state,*E;

	if (nargin<1||mxGetM(argin[0])<6) mexErrMsgTxt("argin error");
	if (nargout>1) mexErrMsgTxt("argout error"); 
	state=mxGetPr(argin[0]);
	argout[0]=mxCreateDoubleMatrix(3,3,mxREAL); E=mxGetPr(argout[0]);
	EcsfToSatf(state,E);
}
