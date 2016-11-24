/*------------------------------------------------------------------------------
% [system] : GpsTools
% [module] : double difference
% [func]   : generate double difference between satellites/stations
% [argin]  : dz  = residuals (O-C) (m)
%            dgds= partial derivatives
%            sig = measurement noise std. deviation (m)
%            ig  = observation index [t,isat,ircv,arcf;...]
%            ircv= station numbers of baselines [ircv1,ircv2;...]
% [argout] : dzs = double diff. of resuduals
%            G   = double diff. of partial derivatives
%            R   = covariance matrix of measurement noise
%            iz  = satellite/station index[isat1,isat2,ircv1,ircv2;...]
% [note]   :
% [version]: $Revision: 7 $ $Date: 06/07/08 1:01 $
%            Copyright(c) 2004-2006 by T.Takasu, all rights reserved
% [history]: 04/06/19  0.1  new
%-----------------------------------------------------------------------------*/
#include "mex.h"

#define		NSATMAX		32		/* max satellite number */
#define		NZMAX		1024	/* max observation data */
#define		NZSMAX		512		/* max double difference */

extern int DblDiff(const double *dz, const double *dgds, const double *sig,
          	       const double *ig, const double *ircv, int nz, int nx, int nr,
          	       double *dzs, double *G, double *R, double *iz);

/* mex interface -------------------------------------------------------------*/
extern void mexFunction(int nargout, mxArray *argout[], int nargin,
					    const mxArray *argin[])
{
	int i,j,nz,nx,nr,nzs;
	double *dz,*dgds,*sig,*ig,*ircv,*_dzs,*_G,*_R,*_iz,*dzs,*G,*R,*iz;

	if (nargin<5||mxGetM(argin[4])<1||mxGetN(argin[4])<2)
		mexErrMsgTxt("argin error");
    nz=mxGetM(argin[0]); if (nz>NZMAX) mexErrMsgTxt("obs data size overflow\n");
    nx=mxGetN(argin[1]);
    nr=mxGetM(argin[4]);
	if (mxGetM(argin[1])!=nz||mxGetM(argin[2])!=nz||mxGetM(argin[3])!=nz)
		mexErrMsgTxt("argin error");
	if (nargout<3) mexErrMsgTxt("argout error");
	if (nz<=0) {
		argout[0]=mxCreateDoubleMatrix(0,1,mxREAL);
		argout[1]=mxCreateDoubleMatrix(0,nx,mxREAL);
		argout[2]=mxCreateDoubleMatrix(0,0,mxREAL);
		if (nargout>=4) argout[3]=mxCreateDoubleMatrix(0,4,mxREAL);
		return;
	}
	dz  =mxGetPr(argin[0]);
	dgds=mxGetPr(argin[1]);
	sig =mxGetPr(argin[2]);
	ig  =mxGetPr(argin[3]);
	ircv=mxGetPr(argin[4]);
	
	_dzs=(double *)malloc(sizeof(double)*NZSMAX);
	_G  =(double *)malloc(sizeof(double)*NZSMAX*nx);
	_R  =(double *)malloc(sizeof(double)*NZSMAX*NZSMAX);
	_iz =(double *)malloc(sizeof(double)*NZSMAX*4);
		
	nzs=DblDiff(dz,dgds,sig,ig,ircv,nz,nx,nr,_dzs,_G,_R,_iz);

	argout[0]=mxCreateDoubleMatrix(nzs,1,mxREAL);   dzs=mxGetPr(argout[0]);
	argout[1]=mxCreateDoubleMatrix(nzs,nx,mxREAL);  G  =mxGetPr(argout[1]);
	argout[2]=mxCreateDoubleMatrix(nzs,nzs,mxREAL); R  =mxGetPr(argout[2]);
	for (i=0;i<nzs;i++) dzs[i]=_dzs[i];
	for (i=0;i<nx; i++) for (j=0;j<nzs;j++) G[nzs*i+j]=_G[NZSMAX*i+j];
	for (i=0;i<nzs;i++) for (j=0;j<nzs;j++) R[nzs*i+j]=_R[NZSMAX*i+j];
	if (nargout>=4) {
		argout[3]=mxCreateDoubleMatrix(nzs,4,mxREAL); iz=mxGetPr(argout[3]);
		for (i=0;i<4;i++) for (j=0;j<nzs;j++) iz[nzs*i+j]=_iz[NZSMAX*i+j];
	}
	free(_dzs); free(_G); free(_R); free(_iz);
}
