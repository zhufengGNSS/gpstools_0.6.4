/*------------------------------------------------------------------------------
% [system] : GpsTools
% [module] : state transition model - receiver clock
% [func]   : state transition model by receiver clock
% [argin]  : t0  = initial time (mjd-utc)
%            t   = time (sec)
%            x0  = initial state(position/velocity) [r0;v0];
%            ircv= staion index (no use)
% [argout] : x   = transit state [r;v]
%            Phi = state transition matrix
% [note]   :
% [version]: $Revision: 4 $ $Date: 06/07/08 1:01 $
%            Copyright(c) 2004-2006 by T.Takasu, all rights reserved
% [history]: 2004/01/05  0.1  new
%-----------------------------------------------------------------------------*/
#include "mex.h"
#include "qtcmn.h"

/* state transition model - receiver clock -----------------------------------*/
extern void State_RcvClock(const double *t0, const double *t, const double *x0,
						   const double *ircv, double *x, double *Phi)
{
	Phi[0]=Phi[3]=1.0;
	Phi[2]=t[0];
	x[0]=x0[0]+t[0]*x0[1];
	x[1]=x0[1];
}
/* mex interface -------------------------------------------------------------*/
extern void mexFunction(int nargout, mxArray *argout[], int nargin,
					    const mxArray *argin[])
{
	double *t0,*t,*x0,*ircv,*x,*Phi;

    if (nargin<4||mxGetN(argin[0])<1||mxGetN(argin[1])<1||mxGetM(argin[2])<2)
        mexErrMsgTxt("argin error");
	if (nargout>2) mexErrMsgTxt("argout error"); 
    
    t0  =mxGetPr(argin[0]);
    t   =mxGetPr(argin[1]);
    x0  =mxGetPr(argin[2]);
    ircv=mxGetPr(argin[3]);
	argout[0]=mxCreateDoubleMatrix(2,1,mxREAL); x  =mxGetPr(argout[0]);
	argout[1]=mxCreateDoubleMatrix(2,2,mxREAL); Phi=mxGetPr(argout[1]);
	
	State_RcvClock(t0,t,x0,ircv,x,Phi);
}
