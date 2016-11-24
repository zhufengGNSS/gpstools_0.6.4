/*------------------------------------------------------------------------------
% [system] : GpsTools
% [module] : ode113 interpolation
% [func]   : ode113 interpolation
% [argin]  : tinterp,t,y,tnew,ynew,klast,phi,psi = refer ntrp113.m
% [argout] : yinterp = refer ntrp113.m
% [note]   : for ode113 speedup (matlab 5.3)
% [version]: $Revision: 3 $ $Date: 06/07/08 1:01 $
%            Copyright(c) 2004-2006 by T.Takasu, all rights reserved
% [history]: 04/01/05  0.1  new
%-----------------------------------------------------------------------------*/
#include "mex.h"

/* ntrp113 -------------------------------------------------------------------*/
extern void ntrp113(const double *tinterp, const double *t, const double *y,
				    const double *tnew, const double *ynew, const double *klast,
				    const double *phi, const double *psi, int np, double *yinterp)
{
	int i,j,ki;
	double hi,w[13],g[13],gamma,eta,term,pg;
	
	hi=tinterp[0]-tnew[0];
	ki=(int)klast[0]+1;
	for (i=0;i<13;i++) {
		w[i]=1.0/(double)(i+1);
		g[i]=0.0;
	}
	g[0]=1.0;
	term=0.0;
	for (j=1;j<ki;j++) {
		gamma=(hi+term)/psi[j-1];
		eta=hi/psi[j-1];
		for (i=0;i<ki+1-j;i++) w[i]=gamma*w[i]-eta*w[i+1];
		g[j]=w[0];
		term=psi[j-1];
	}
	for (j=0;j<np;j++) {
		for (i=0,pg=0.0;i<ki;i++) pg+=phi[j+i*np]*g[i];
		yinterp[j]=ynew[j]+hi*pg;
	}
}
/* mex interface ------------------------------------------------------------*/
extern void mexFunction(int nargout, mxArray *argout[], int nargin,
					    const mxArray *argin[])
{
	int np;
	double *tinterp,*t,*y,*tnew,*ynew,*klast,*phi,*psi,*yinterp;

	if (nargin<8) mexErrMsgTxt("argin error");
	if (nargout>1) mexErrMsgTxt("argout error");
	np=mxGetM(argin[6]);
	tinterp=mxGetPr(argin[0]);
	t      =mxGetPr(argin[1]);
	y      =mxGetPr(argin[2]);
	tnew   =mxGetPr(argin[3]);
	ynew   =mxGetPr(argin[4]);
	klast  =mxGetPr(argin[5]);
	phi    =mxGetPr(argin[6]);
	psi    =mxGetPr(argin[7]);
/*
	mexPrintf("tinterp(%d,%d)=%f\n",mxGetM(argin[0]),mxGetN(argin[0]),*tinterp);
	mexPrintf("t(%d,%d)\n",mxGetM(argin[1]),mxGetN(argin[1]));
	mexPrintf("y(%d,%d)\n",mxGetM(argin[2]),mxGetN(argin[2]));
	mexPrintf("tnew(%d,%d)=%f\n",mxGetM(argin[3]),mxGetN(argin[3]),*tnew);
	mexPrintf("ynew(%d,%d)\n",mxGetM(argin[4]),mxGetN(argin[4]));
	mexPrintf("klast(%d,%d)=%f\n",mxGetM(argin[5]),mxGetN(argin[5]),*klast);
	mexPrintf("phi(%d,%d)\n",mxGetM(argin[6]),mxGetN(argin[6]));
	mexPrintf("psi(%d,%d)\n",mxGetM(argin[7]),mxGetN(argin[7]));
	mexErrMsgTxt("");
*/
	argout[0]=mxCreateDoubleMatrix(np,1,mxREAL); yinterp=mxGetPr(argout[0]);
	ntrp113(tinterp,t,y,tnew,ynew,klast,phi,psi,np,yinterp);
}
