/*------------------------------------------------------------------------------
% [system] : GpsTools
% [module] : relativity correction
% [func]   : relativity correction
% [argin]  : rsat = satellite postion(m) (eci)
%            vsat = satellite velocity(m/sec) (eci)
%            rrcv = station position(m) (ecef)
%           (opts)= option flag (1=shapiro-delay correction) (default:1)
% [argout] : rels = relativity correction(m)
% [note]   :
% [version]: $Revision: 16 $ $Date: 2008-12-12 15:49:30 +0900 (金, 12 12 2008) $
%            Copyright(c) 2004-2006 by T.Takasu, all rights reserved
% [history]: 04/10/19   0.1  new
%            08/11/25   0.2  add argin opts (gt_0.6.4)
%-----------------------------------------------------------------------------*/
#include "mex.h"

extern double RelCorr(const double *rsat, const double *vsat,
					  const double *rrcv, int opts);

/* mex inerface --------------------------------------------------------------*/
extern void mexFunction(int nargout, mxArray *argout[], int nargin,
					    const mxArray *argin[])
{
	double *rsat,*vsat,*rrcv,*rels;
	int opts=1;
	
	if (nargin<3||mxGetM(argin[0])<3||mxGetM(argin[1])<3||mxGetM(argin[2])<3)
		mexErrMsgTxt("argin error");
	
	rsat=mxGetPr(argin[0]);
	vsat=mxGetPr(argin[1]);
	rrcv=mxGetPr(argin[2]);
	if (nargin>3) {
		opts=(int)*mxGetPr(argin[3]);
	}
	argout[0]=mxCreateDoubleMatrix(1,1,mxREAL); rels=mxGetPr(argout[0]);
	
	*rels=RelCorr(rsat,vsat,rrcv,opts);
}
