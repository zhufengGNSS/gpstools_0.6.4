/*------------------------------------------------------------------------------
% [system] : GpsTools
% [module] : spheric harmonic functions
% [func]   : spheric harmonic functions by azimuth/elevation angle
% [argin]  : az,el = azimuth/elevation angle (rad)
%            nmax  = max degrees of spheric harmonic functions
% [argout] : fc,fs = spheric harmonic functions
%                fc(n+1,m+1) = Pnm(-cos(2*el))*cos(m*az)
%                fs(n+1,m+1) = Pnm(-cos(2*el))*sin(m*az)
%                (Pnm = normalized Legendre polynomial)
% [note]   :
% [version]: $Revision: 4 $ $Date: 06/07/08 1:01 $
%            Copyright(c) 2004-2006 by T.Takasu, all rights reserved
% [history]: 04/10/18  0.1  new
%-----------------------------------------------------------------------------*/
#include "mex.h"

extern void SphFunc(double az, double el, int nmax, double *fc, double *fs);

/* mex interface -------------------------------------------------------------*/
extern void mexFunction(int nargout, mxArray *argout[], int nargin,
					    const mxArray *argin[])
{
	double az,el,*fc,*fs;
	int nmax;
	
	if (nargin<3||mxGetM(argin[0])<1||mxGetM(argin[1])<1||mxGetM(argin[2])<1)
		mexErrMsgTxt("argin error");
	
	az=*mxGetPr(argin[0]);
	el=*mxGetPr(argin[1]);
	nmax=(int)*mxGetPr(argin[2]);
	argout[0]=mxCreateDoubleMatrix(nmax+1,nmax+1,mxREAL); fc=mxGetPr(argout[0]);
	argout[1]=mxCreateDoubleMatrix(nmax+1,nmax+1,mxREAL); fs=mxGetPr(argout[1]);
	
	SphFunc(az,el,nmax,fc,fs);
}
