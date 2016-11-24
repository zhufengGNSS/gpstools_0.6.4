/*------------------------------------------------------------------------------
% [system] : GpsTools
% [module] : state transition model - two-body problem
% [func]   : state transition model by two-body problem (analytical)
% [argin]  : t0  = initial time (mjd-utc)
%            t   = time (sec)
%            x0  = initial state(position/velocity) [r0;v0];
%            isat= satellite index (no use)
% [argout] : x   = transit state [r;v]
%            Phi = state transition matrix
% [note]   :
% [version]: $Revision: 12 $ $Date: 2008-11-25 10:02:15 +0900 (火, 25 11 2008) $
%            Copyright(c) 2004-2006 by T.Takasu, all rights reserved
% [history]: 2004/01/05  0.1  new
%-----------------------------------------------------------------------------*/
#include "mex.h"
#include "qtcmn.h"

extern void StateToEle(const double *state, double *ele);
extern void EleToState(const double *ele, double *state, double *dsde);

/* state transition model - two-body problem (analytical) --------------------*/
extern void State_Kepler(const double *t0, const double *t, const double *x0,
						 const double *isat, double *x, double *Phi)
{
	int i;
	double a3,a5,xx[6],ele[6],ds0de0[36],dsde[36],ids0de0[36],dPhi[36];
	
	StateToEle(x0,ele);
	a3=ele[0]*ele[0]*ele[0]; a5=a3*ele[0]*ele[0];
	EleToState(ele,xx,ds0de0);
	ele[5]+=sqrt(GME/a3)*t[0]*180.0/M_PI;
	EleToState(ele,x,dsde);
	for (i=0;i<36;i+=7) Phi[i]=1.0;
	Phi[5]=-1.5*sqrt(GME/a5)*t[0];
	MatInv(ds0de0,6,ids0de0);
	MatMul(dsde,Phi,6,6,6,dPhi);
	MatMul(dPhi,ids0de0,6,6,6,Phi);
}
/* mex interface -------------------------------------------------------------*/
extern void mexFunction(int nargout, mxArray *argout[], int nargin,
					    const mxArray *argin[])
{
	int i,j,nx;
	double *t0,*t,*x0,*isat,*x,*Phi,Phix[36];

    if (nargin<4||mxGetN(argin[0])<1||mxGetN(argin[1])<1||mxGetM(argin[2])<6)
        mexErrMsgTxt("argin error");
	if (nargout!=2) mexErrMsgTxt("argout error"); 
    
    t0  =mxGetPr(argin[0]);
    t   =mxGetPr(argin[1]);
    x0  =mxGetPr(argin[2]); nx=mxGetM(argin[2]);
    isat=mxGetPr(argin[3]);
	argout[0]=mxCreateDoubleMatrix(nx,1,mxREAL); x=mxGetPr(argout[0]);
	argout[1]=mxCreateDoubleMatrix(nx,nx,mxREAL); Phi=mxGetPr(argout[1]);
	
	State_Kepler(t0,t,x0,isat,x,Phix);
	for (i=0;i<6;i++) for (j=0;j<6;j++) Phi[i*nx+j]=Phix[i*6+j];
	for (i=6;i<nx;i++) {x[i]=x0[i]; Phi[i*nx+i]=1.0;}
}
