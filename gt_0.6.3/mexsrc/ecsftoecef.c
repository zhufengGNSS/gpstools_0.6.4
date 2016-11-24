/*------------------------------------------------------------------------------
% [system] : GpsTools
% [module] : eci to ecef transformation matrix
% [func]   : calculate eci to ecef transformation matrix
% [argin]  : tm = date/time (mjd-utc)
%            (erp_value)= earth rotation parameters (default:[0,0,0,0,0])
%                erpvalue(1:2) = pole offset xp,yp(rad)
%                erpvalue(3)   = ut1-utc(sec)
%                erpvalue(4:5) = nutation offset dpsi,deps(rad)
%            (utc_tai)  = utc-tai(sec) (default:-32)
%            (model_nut)= nutation model
%                ''         : iau1976/iau1980+dpsi,deps
%                'iers1996' : iers1996
% [argout] : U    = eci to ecef transformation matrix(3x3)
%            P,N  = precession,nutation matrix(3x3)
%            gmst = Greenwich mean sidereal time(rad)
%            (dUdxp,dUdyp,dUddt) = dU/dxp,dU/dyp,dU/d(UT1-UTC)(1/rad,1/rad,1/sec)
% [note]   : iau1976 precession,iau1980 nutation model
%            eci  : icrf (J2000.0 mean equator coordinate)
%            ecef : itrf
%            inv(U)=U',inv(P)=P',inv(N)=N'
% [version]: $Revision: 6 $ $Date: 06/07/08 1:01 $
%            Copyright(c) 2004-2006 by T.Takasu, all rights reserved
% [history]: 03/12/31  0.1  new
%-----------------------------------------------------------------------------*/
#include "mex.h"

extern void EcsfToEcef(const double *tm, const double *erp_value,
                       const double *utc_tai,const char *model_nut, double *U,
                       double *P, double *N, double *gmst, double *dUdxp,
                       double *dUdyp, double *dUddt, double *fargs);

/* mex interface -------------------------------------------------------------*/
extern void mexFunction(int nargout, mxArray *argout[], int nargin,
					    const mxArray *argin[])
{
	double erp_value_v[]={0.0,0.0,0.0,0.0,0.0},utc_tai_v[]={-32.0};
	double *erp_value=erp_value_v,*utc_tai=utc_tai_v,Pv[9],Nv[9],gmstv[1];
	double *tm,*U,*P=Pv,*N=Nv,*gmst=gmstv,*dUdxp=NULL,*dUdyp=NULL,*dUddt=NULL;
	char model_nut[32]="";

    if (nargin<1||mxGetN(argin[0])<1) mexErrMsgTxt("argin error");
	if (nargout>7) mexErrMsgTxt("argout error"); 
    
    tm=mxGetPr(argin[0]);
	if (nargin>1&&mxGetN(argin[1])>=5) erp_value=mxGetPr(argin[1]);
	if (nargin>2&&mxGetN(argin[2])>=1) utc_tai=mxGetPr(argin[2]);
	if (nargin>3&&mxIsChar(argin[3])) mxGetString(argin[3],model_nut,32);
	
	argout[0]=mxCreateDoubleMatrix(3,3,mxREAL); U=mxGetPr(argout[0]);
	if (nargout>1) {argout[1]=mxCreateDoubleMatrix(3,3,mxREAL); P=mxGetPr(argout[1]);}
	if (nargout>2) {argout[2]=mxCreateDoubleMatrix(3,3,mxREAL); N=mxGetPr(argout[2]);}
	if (nargout>3) {argout[3]=mxCreateDoubleMatrix(1,1,mxREAL); gmst=mxGetPr(argout[3]);}
	if (nargout>4) {argout[4]=mxCreateDoubleMatrix(3,3,mxREAL); dUdxp=mxGetPr(argout[4]);}
	if (nargout>5) {argout[5]=mxCreateDoubleMatrix(3,3,mxREAL); dUdyp=mxGetPr(argout[5]);}
	if (nargout>6) {argout[6]=mxCreateDoubleMatrix(3,3,mxREAL); dUddt=mxGetPr(argout[6]);}
    
    EcsfToEcef(tm,erp_value,utc_tai,model_nut,U,P,N,gmst,dUdxp,dUdyp,dUddt,NULL);
}
