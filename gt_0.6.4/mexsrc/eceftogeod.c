/*------------------------------------------------------------------------------
% [system] : GpsTools
% [module] : ecef to geodetic position
% [func]   : convert ecef postion to geodetic (latitude,longitude,height)
% [argin]  : epos = ecef position (m) [x;y;z]
% [argout] : gpos = latitude/longitude/height(deg,m) [lat,lon,h] 
%           (E)   = ecef to local tangental coordinate transformation matrix
% [note]   :
% [version]: $Revision: 12 $ $Date: 2008-11-25 10:02:15 +0900 (火, 25 11 2008) $
%            Copyright(c) 2004-2006 by T.Takasu, all rights reserved
% [history]: 04/05/31   0.1  new
%            05/12/12   0.2  add argout E
%-----------------------------------------------------------------------------*/
#include "mex.h"

extern void EcefToGeod(const double *epos, double *gpos, double *E);

/* mex interface -------------------------------------------------------------*/
extern void mexFunction(int nargout, mxArray *argout[], int nargin,
					    const mxArray *argin[])
{
	double *epos,*gpos,*E=NULL;

	if (nargin<1||mxGetM(argin[0])<3) mexErrMsgTxt("argin error");
	if (nargout>2) mexErrMsgTxt("argout error"); 
	epos=mxGetPr(argin[0]);
	argout[0]=mxCreateDoubleMatrix(1,3,mxREAL); gpos=mxGetPr(argout[0]);
	if (nargout>1) {
		argout[1]=mxCreateDoubleMatrix(3,3,mxREAL); E=mxGetPr(argout[1]);
	}
	EcefToGeod(epos,gpos,E);
}
