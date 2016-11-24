/*------------------------------------------------------------------------------
% [system] : GpsTools
% [module] : tropospheric mapping function - NMF
% [func]   : calculate tropospheric mapping function by Niell mapping function
% [argin]  : t    = date/time (mjd)
%            azel = azimath/elevation angle(rad) [az,el]
%            gpos = latitude/longitude/height(deg,m) [lat,lon,h]
%            (metprm) = meteological parameters(no use)
% [argout] : mapfd = dry mapping function
%            mapfw = wet mapping function
% [note]   : reference :
%            A.E.Niell, Global mapping functions for the atmosphere delay at
%            radio wavelengths, 1996
% [version]: $Revision: 12 $ $Date: 2008-11-25 10:02:15 +0900 (火, 25 11 2008) $
%            Copyright(c) 2004-2006 by T.Takasu, all rights reserved
% [history]: 04/05/31   0.1  new
%-----------------------------------------------------------------------------*/
#include "mex.h"

extern void mapf_nmf(double *t, const double *azel, const double *gpos,
					 double *mapfd, double *mapfw);

/* mex interface -------------------------------------------------------------*/
extern void mexFunction(int nargout, mxArray *argout[], int nargin,
					    const mxArray *argin[])
{
	double *t,*azel,*gpos,*mapfd,*mapfw;
	
	if (nargin<3||mxGetM(argin[0])<1||mxGetN(argin[1])<2||mxGetN(argin[2])<3)
		mexErrMsgTxt("argin error");
	
	t   =mxGetPr(argin[0]);
	azel=mxGetPr(argin[1]);
	gpos=mxGetPr(argin[2]);
	argout[0]=mxCreateDoubleMatrix(1,1,mxREAL); mapfd=mxGetPr(argout[0]);
	argout[1]=mxCreateDoubleMatrix(1,1,mxREAL); mapfw=mxGetPr(argout[1]);
	
	mapf_nmf(t,azel,gpos,mapfd,mapfw);
}
