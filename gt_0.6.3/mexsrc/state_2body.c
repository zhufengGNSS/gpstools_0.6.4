/*------------------------------------------------------------------------------
% [system] : GpsTools
% [module] : state transition model - two-body problem
% [func]   : state transition model by two-body problem (numerical integration)
% [argin]  : t0  = initial time (mjd-utc)
%            t   = time (sec)
%            x0  = initial state(position/velocity) [r0;v0];
%            isat= satellite index (no use)
% [argout] : x   = transit state [r;v]
%            Phi = state transition matrix
% [note]   :
% [version]: $Revision: 4 $ $Date: 06/07/08 1:01 $
%            Copyright(c) 2004-2006 by T.Takasu, all rights reserved
% [history]: 2004/01/05  0.1  new
%-----------------------------------------------------------------------------*/
#include "mex.h"
#include "qtcmn.h"

#define NX			(6+6*6)				/* number of states */
#define TSTEP		1.0					/* integration step time (sec) */

/* ordinary differential equation --------------------------------------------*/
static void deq(double t, const double *x, double *xdot)
{
	int i;
	double r2=DOT(x,x),r3=r2*sqrt(r2),dadr[9];

	for (i=0;i<3;i++) {
		xdot[i]=x[i+3];
		xdot[i+3]=-GME*x[i]/r3;
	}
	dadr[0]=3.0*x[0]*x[0]-r2;
	dadr[4]=3.0*x[1]*x[1]-r2;
	dadr[8]=3.0*x[2]*x[2]-r2;
	dadr[1]=dadr[3]=3.0*x[0]*x[1];
	dadr[5]=dadr[7]=3.0*x[1]*x[2];
	dadr[2]=dadr[6]=3.0*x[2]*x[0];
	for (i=0;i<9;i++) {
		dadr[i]*=GME/(r2*r3);
		xdot[i+6]=x[i+15];		/* dƒ³rr=ƒ³vr */
		xdot[i+24]=x[i+33];		/* dƒ³rv=ƒ³vv */
	}
	MM(dadr,x+6,xdot+15);		/* dƒ³vr=dadr*ƒ³rr */
	MM(dadr,x+24,xdot+33);		/* dƒ³vv=dadr*ƒ³rv */
}
/* numerical integration -----------------------------------------------------*/
static void RK4(double t, double *x)
{
	int i;
	double tt,h,k1[NX],k2[NX],k3[NX],k4[NX],w[NX];

	for (tt=0.0;tt<t-1E-12;tt=tt+h) {
	    h=MIN(t-tt,TSTEP);
	    deq(tt,x,k1);	  for (i=0;i<NX;i++) w[i]=x[i]+k1[i]*h/2.0;
	    deq(tt+h/2,w,k2); for (i=0;i<NX;i++) w[i]=x[i]+k2[i]*h/2.0;
	    deq(tt+h/2,w,k3); for (i=0;i<NX;i++) w[i]=x[i]+k3[i]*h;
	    deq(tt+h,w,k4);
		for (i=0;i<NX;i++) x[i]=x[i]+(k1[i]+2.0*k2[i]+2.0*k3[i]+k4[i])*h/6.0;
	}
}
/* state transition model - two-body problem ---------------------------------*/
extern void State_2Body(const double *t0, const double *t, const double *x0,
						const double *isat, double *x, double *Phi)
{
	static const int index[]={
		0,1,2,9,10,11, 3,4,5,12,13,14, 6,7,8,15,16,17,
		18,19,20,27,28,29, 21,22,23,30,31,32, 24,25,26,33,34,35
	};
	int i;
	double xx[NX];

	for (i=0;i<6;i++) xx[i]=x0[i];
	EYE(xx+6); ZERO(xx+15); ZERO(xx+24); EYE(xx+33);
	
	RK4(t[0],xx);
	
	for (i=0;i<6;i++) x[i]=xx[i];
	for (i=0;i<36;i++) Phi[i]=xx[index[i]+6];
}
/* mex interface -------------------------------------------------------------*/
extern void mexFunction(int nargout, mxArray *argout[], int nargin,
					    const mxArray *argin[])
{
	int i,j,nx;
	double *t0,*t,*x0,*isat,*x,*Phi,Phix[36];

    if (nargin<4||mxGetN(argin[0])<1||mxGetN(argin[1])<1||mxGetM(argin[2])<6)
        mexErrMsgTxt("argin error");
	if (nargout>2) mexErrMsgTxt("argout error"); 
    
    t0  =mxGetPr(argin[0]);
    t   =mxGetPr(argin[1]);
    x0  =mxGetPr(argin[2]); nx=mxGetM(argin[2]);
    isat=mxGetPr(argin[3]);
	argout[0]=mxCreateDoubleMatrix(nx,1,mxREAL);  x  =mxGetPr(argout[0]);
	argout[1]=mxCreateDoubleMatrix(nx,nx,mxREAL); Phi=mxGetPr(argout[1]);
	
	State_2Body(t0,t,x0,isat,x,Phix);
	for (i=0;i<6;i++) for (j=0;j<6;j++) Phi[i*nx+j]=Phix[i*6+j];
	for (i=6;i<nx;i++) {x[i]=x0[i]; Phi[i*nx+i]=1.0;}
}
