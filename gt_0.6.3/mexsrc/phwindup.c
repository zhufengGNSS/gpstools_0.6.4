/*------------------------------------------------------------------------------
% [system] : GpsTools
% [module] : phase windup effect correction
% [func]   : calculate phase windup effect correction
% [argin]  : rsat = satellite position(m) (eci)
%            rsun = sun position(m) (eci)
%            posr = station position (m) (ecef)
%            U    = eci to ecef transformation matrix
%            phwinp = previous correction(rad)
%           (rmode) = rcv attitude mode (0:fixed station,1:leo satellite)
% [argout] : phwin  = phase windup phase correction(rad)
% [note]   :
% [version]: $Revision: 7 $ $Date: 06/07/08 1:01 $
%            Copyright(c) 2004-2006 by T.Takasu, all rights reserved
% [history]: 04/06/18   0.1  new
%            04/08/12   0.2  add argin rmode
%-----------------------------------------------------------------------------*/
#include "mex.h"

extern double phwindup(const double *rsat, const double *rsun,
					   const double *posr, const double *U, const double *phwinp,
					   const double *rmode);

/* mex interface -------------------------------------------------------------*/
extern void mexFunction(int nargout, mxArray *argout[], int nargin,
					    const mxArray *argin[])
{
	double *rsat,*rsun,*posr,*U,*phwinp,*phwin,*rmode,rm=0.0;
	
	if (nargin<5||mxGetM(argin[0])<3||mxGetM(argin[1])<3||mxGetM(argin[2])<3||
	    mxGetM(argin[3])!=3||mxGetN(argin[3])!=3||mxGetN(argin[4])<1)
		mexErrMsgTxt("argin error");
	
	rsat  =mxGetPr(argin[0]);
	rsun  =mxGetPr(argin[1]);
	posr  =mxGetPr(argin[2]);
	U     =mxGetPr(argin[3]);
	phwinp=mxGetPr(argin[4]);
	if (nargin>=6) rmode=mxGetPr(argin[5]); else rmode=&rm;
	argout[0]=mxCreateDoubleMatrix(1,1,mxREAL); phwin=mxGetPr(argout[0]);
	
	*phwin=phwindup(rsat,rsun,posr,U,phwinp,rmode);
}
