/*------------------------------------------------------------------------------
% [system] : GpsTools
% [module] : state transition model - j2 perturbation
% [func]   : state transition model by j2 perturbation model
% [argin]  : t0  = initial time (mjd-utc)
%            t   = time (sec)
%            x0  = initial state(position/velocity) [r0;v0];
%            isat= satellite index (no use)
% [argout] : x   = transit state [r;v]
%            Phi = state transition matrix
% [note]   :
% [version]: $Revision: 5 $ $Date: 06/07/08 1:01 $
%            Copyright(c) 2004-2006 by T.Takasu, all rights reserved
% [history]: 2004/01/05  0.1  new
%-----------------------------------------------------------------------------*/
#include "mex.h"
#include "qtcmn.h"

#define NX			(6+6*6)				/* number of states */
#define TSTEP		1.0					/* integration step time (sec) */

static double Ugeoj2[9];

extern void EcsfToEcef(const double *tm, const double *erp_value,
                       const double *utc_tai,const char *model_nut, double *U,
                       double *P, double *N, double *gmst, double *dUdxp,
                       double *dUdyp, double *dUddt, double *fargs);

/* geopotential acceleration and partial derivatives with j2 -----------------*/
static void Accel(const double *r, double *a, double *dadr)
{
    static const double GMe=3.986004415E+14;
    static const double Re=6378136.3;
    static const double J2=0.001082636;
    double k=1.5*J2*Re*Re;
    double r2,r3,p1,p2,daxdx,daydy,daxdy,daxdz,daydz;

    r2=DOT(r,r); r3=r2*sqrt(r2);
    p1=1.0-5.0*r[2]*r[2]/r2;
    p2=3.0+5.0*k/r2-35.0*k*r[2]*r[2]/(r2*r2);
    a[0]=-GMe/r3*(r[0]+k/r2*p1*r[0]);
    a[1]=-GMe/r3*(r[1]+k/r2*p1*r[1]);
    a[2]=-GMe/r3*(r[2]+k/r2*(p1+2.0)*r[2]);
    daxdx=GMe/(r2*r3)*(r[0]*r[0]*p2-k*p1-r2);
    daydy=GMe/(r2*r3)*(r[1]*r[1]*p2-k*p1-r2);
    daxdy=GMe/(r2*r3)*r[0]*r[1]*p2;
    daxdz=GMe/(r2*r3)*r[0]*r[2]*(p2+10.0*k/r2);
    daydz=GMe/(r2*r3)*r[1]*r[2]*(p2+10.0*k/r2);
    dadr[0]=daxdx;
    dadr[1]=dadr[3]=daxdy;
    dadr[2]=dadr[6]=daxdz;
    dadr[4]=daydy;
    dadr[5]=dadr[7]=daydz;
    dadr[8]=-daxdx-daydy;
}
/* ordinary differential equation --------------------------------------------*/
static void deq(double t, const double *x, double *xdot)
{
	int i;
	double Ut[9],r[3],a[3],Ut_a[3],dadr[9],Ut_dadr[9];
    
    Tr(Ugeoj2,Ut);
    Mv(Ugeoj2,x,r);
    Accel(r,a,dadr);
    Mv(Ut,a,Ut_a);
    MM(Ut,dadr,Ut_dadr);
    MM(Ut_dadr,Ugeoj2,dadr);

	for (i=0;i<3;i++) {
		xdot[i]=x[i+3];
		xdot[i+3]=Ut_a[i];
	}
	for (i=0;i<9;i++) {
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
/* state transition model - j2 perturbation ----------------------------------*/
extern void State_GeoJ2(const double *t0, const double *t, const double *x0,
						const double *isat, double *x, double *Phi)
{
	static const int index[]={
		0,1,2,9,10,11, 3,4,5,12,13,14, 6,7,8,15,16,17,
		18,19,20,27,28,29, 21,22,23,30,31,32, 24,25,26,33,34,35
	};
	int i;
	double xx[NX],P[9],N[9],gmst[1];
    
    EcsfToEcef(t0,NULL,NULL,"",Ugeoj2,P,N,gmst,NULL,NULL,NULL,NULL);

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
	
	State_GeoJ2(t0,t,x0,isat,x,Phix);
	for (i=0;i<6;i++) for (j=0;j<6;j++) Phi[i*nx+j]=Phix[i*6+j];
	for (i=6;i<nx;i++) {x[i]=x0[i]; Phi[i*nx+i]=1.0;}
}
