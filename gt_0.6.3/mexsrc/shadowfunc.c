/*-------------------------------------------------------------------------------
% [system] : GpsTools
% [module] : shadow function
% [func]   : calculate shadow factor in eclipse
% [argin]  : rpos  = satellite position (m)
%            rsun  = sun position (m)
%            rmoon = moon position (m)
% [argout] : fact  = shadow factor (0:in eclipse,1:not in eclipse)
%           (model)= shadow model
%                    'cylind'  : cylindric model
%                    'penumbra': penumbra/umbra model
% [note]   :
% [version]: $Revision: 4 $ $Date: 06/07/08 1:01 $
%            Copyright(c) 2004-2006 by T.Takasu, all rights reserved
% [history]: 04/09/28  0.1  new
%            05/04/01  0.2  add penumbra/umbra model
%-----------------------------------------------------------------------------*/
#include "mex.h"

extern double ShadowFunc(const double *rpos, const double *rsun,
					     const double *rmoon, const char *model);

/* mex interface -------------------------------------------------------------*/
extern void mexFunction(int nargout, mxArray *argout[], int nargin,
					    const mxArray *argin[])
{
	char model[32]="";
	double *rpos,*rsun,*rmoon,*fact;

	if (nargin<3||mxGetM(argin[0])<3||mxGetM(argin[1])<3||mxGetM(argin[2])<3)
		mexErrMsgTxt("argin error");
	if (nargout<1) mexErrMsgTxt("argout error"); 
	
	rpos =mxGetPr(argin[0]);
	rsun =mxGetPr(argin[1]);
	rmoon=mxGetPr(argin[2]);
	if (nargin>=4&&mxGetString(argin[3],model,sizeof(model))!=0)
		 mexErrMsgTxt("argin error");
	
	argout[0]=mxCreateDoubleMatrix(1,1,mxREAL); fact=mxGetPr(argout[0]);
	
	*fact=ShadowFunc(rpos,rsun,rmoon,model);
}
