/*------------------------------------------------------------------------------
% [system] : GpsTools
% [module] : geometric distance
% [func]   : geometric distance
% [argin]  : td,t = time mjd,sec
%            nav  = navigation message
%            rr   = station position (m)
% [argout] : r = geometric distance
%            dr= partial derivative of r by rr
% [note]   :
% [version]: $Revision: 4 $ $Date: 06/07/08 1:01 $
%            Copyright(c) 2004-2006 by T.Takasu, all rights reserved
% [history]: 06/01/11  0.1  new
%-----------------------------------------------------------------------------*/
#include "mex.h"
#include "qtcmn.h"

#define OMGE	7.2921151467E-5			/* earth rotation rate (rad/sec) */

extern void NavToState(double td, double ts, const double *nav, int nnav,
					   char type, int opt, double *pos, double *dts, double *vel);

/* geometric distance --------------------------------------------------------*/
static void geodist(const double *td, const double *t, const double *nav,
					int nnav, const double *rr, double *r, double *dr)
{
	int i;
	double tau=0.0,rs[3],rss[3],R[9],rr_rs[3];
	
	while (1) {
		NavToState(*td,*t-tau,nav,nnav,'N',0,rs,NULL,NULL);
		Rz(OMGE*tau,R);
		Mv(R,rs,rss);
		for (i=0;i<3;i++) rr_rs[i]=rr[i]-rss[i];
		*r=NORM(rr_rs);
		if (ABS(*r-tau*M_C)<1E-4) break;
		tau=*r/M_C;
	}
	for (i=0;i<3;i++) dr[i]=rr_rs[i]/(*r);
}
/* mex interface -------------------------------------------------------------*/
extern void mexFunction(int nargout, mxArray *argout[], int nargin,
					    const mxArray *argin[])
{
	int nnav;
	double *td,*t,*nav,*rr,*r,*dr;
	
	if (nargin!=4) mexErrMsgTxt("argin error");
	if (nargout!=2) mexErrMsgTxt("argout error"); 
	
	if (mxGetM(argin[0])<1||mxGetM(argin[1])<1||mxGetN(argin[2])<31||
		mxGetM(argin[3])<3)
		mexErrMsgTxt("argin error");
	
	td=mxGetPr(argin[0]);
	t=mxGetPr(argin[1]);
	nav=mxGetPr(argin[2]); nnav=mxGetM(argin[2]);
	rr=mxGetPr(argin[3]);
	argout[0]=mxCreateDoubleMatrix(1,1,mxREAL); r=mxGetPr(argout[0]);
	argout[1]=mxCreateDoubleMatrix(1,3,mxREAL); dr=mxGetPr(argout[1]);
	geodist(td,t,nav,nnav,rr,r,dr);
}
