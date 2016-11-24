/*------------------------------------------------------------------------------
% [system] : GpsTools
% [module] : relativity correction
% [func]   : relativity correction
% [argin]  : rsat = satellite postion(m) (eci)
%            vsat = satellite velocity(m/sec) (eci)
%            rrcv = station position(m) (ecef)
% [argout] : rels = relativity correction(m)
% [note]   :
% [version]: $Revision: 4 $ $Date: 06/07/08 1:01 $
%            Copyright(c) 2004-2006 by T.Takasu, all rights reserved
% [history]: 04/10/19   0.1  new
%-----------------------------------------------------------------------------*/
#include "mex.h"

extern double RelCorr(const double *rsat, const double *vsat,
					  const double *rrcv);

/* mex inerface --------------------------------------------------------------*/
extern void mexFunction(int nargout, mxArray *argout[], int nargin,
					    const mxArray *argin[])
{
	double *rsat,*vsat,*rrcv,*rels;
	
	if (nargin<3||mxGetM(argin[0])<3||mxGetM(argin[1])<3||mxGetM(argin[2])<3)
		mexErrMsgTxt("argin error");
	
	rsat=mxGetPr(argin[0]);
	vsat=mxGetPr(argin[1]);
	rrcv=mxGetPr(argin[2]);
	argout[0]=mxCreateDoubleMatrix(1,1,mxREAL); rels=mxGetPr(argout[0]);
	
	*rels=RelCorr(rsat,vsat,rrcv);
}
