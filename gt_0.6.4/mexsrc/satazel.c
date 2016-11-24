/*-------------------------------------------------------------------------------
% [system] : GpsTools
% [module] : satellite azimuth/elevation angle
% [func]   : calculate satellite azimuth/elevation angle
% [argin]  : spos = satellite postion [posx;posy;posz] (m) (ecef)
%            rpos = station position [posx;posy;posz] (m) (ecef)
% [argout] : azel = azimuth/elevation angle [az,el] (rad)
% [note]   :
% [version]: $Revision: 12 $ $Date: 2008-11-25 10:02:15 +0900 (火, 25 11 2008) $
%            Copyright(c) 2004-2006 by T.Takasu, all rights reserved
% [history]: 04/05/31   0.1  new
%-----------------------------------------------------------------------------*/
#include "mex.h"

extern void SatAzEl(const double *spos, const double *rpos, double *azel);

/* mex interface -------------------------------------------------------------*/
extern void mexFunction(int nargout, mxArray *argout[], int nargin,
					    const mxArray *argin[])
{
	double *spos,*rpos,*azel;

	if (nargin<2||mxGetM(argin[0])<3||mxGetM(argin[1])<3)
		mexErrMsgTxt("argin error");
	
	spos=mxGetPr(argin[0]);
	rpos=mxGetPr(argin[1]);
	
	argout[0]=mxCreateDoubleMatrix(1,2,mxREAL); azel=mxGetPr(argout[0]);
	
	SatAzEl(spos,rpos,azel);
}
