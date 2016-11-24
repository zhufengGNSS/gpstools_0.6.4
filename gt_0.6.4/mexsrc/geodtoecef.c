/*------------------------------------------------------------------------------
% [system] : GpsTools
% [module] : geodetic to ecef position
% [func]   : convert geodetic postion to ecef
% [argin]  : gpos = latitude/longitude/height(deg,m) [lat,lon,h] 
% [argout] : epos = ecef position (m) [x;y;z]
%            E    = ecef to local tangental coordinate transformation matrix
% [note]   : WGS84 ellipsoide
% [version]: $Revision: 12 $ $Date: 2008-11-25 10:02:15 +0900 (火, 25 11 2008) $
%            Copyright(c) 2004-2006 by T.Takasu, all rights reserved
% [history]: 04/05/31   0.1  new
%-----------------------------------------------------------------------------*/
#include "mex.h"

extern void GeodToEcef(const double *gpos, double *epos, double *E);

/* mex interface -------------------------------------------------------------*/
extern void mexFunction(int nargout, mxArray *argout[], int nargin,
					    const mxArray *argin[])
{
	double *gpos,*epos,*E=NULL;

    if (nargin<1||mxGetN(argin[0])<3) mexErrMsgTxt("argin error");
	if (nargout>2) mexErrMsgTxt("argout error"); 
    
    gpos=mxGetPr(argin[0]);
	argout[0]=mxCreateDoubleMatrix(3,1,mxREAL); epos=mxGetPr(argout[0]);
	if (nargout>1) {argout[1]=mxCreateDoubleMatrix(3,3,mxREAL); E=mxGetPr(argout[1]);}
	
	GeodToEcef(gpos,epos,E);
}
