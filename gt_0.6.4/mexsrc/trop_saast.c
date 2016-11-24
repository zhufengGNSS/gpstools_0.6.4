/*------------------------------------------------------------------------------
% [system] : GpsTools
% [module] : tropospheric model - Saastamoinen model
% [func]   : calculate tropospheric delay by Saastamoinen model
% [argin]  : t    = date/time (mjd-utc)(no use)
%            azel = azimuth/elevation angle(rad) [az,el]
%            gpos = latitude/longitude/height(deg,m) [lat,lon,h]
%           (met_prm) = meteorological parameters
%              met_prm(1) = pressure(hPa)
%              met_prm(2) = temperture(C)
%              met_prm(3) = relative humidity(%)
%              (default:standard atmosphere,humidity:50%)
% [argout] : tropd = total/dry tropospheric delay (m)
%           (tropw)= wet tropospheric delay (m)
% [note]   : Reference: H.Lichtenegger GPS Theory and Practice 6.3
% [version]: $Revision: 12 $ $Date: 2008-11-25 10:02:15 +0900 (火, 25 11 2008) $
%            Copyright(c) 2004-2006 by T.Takasu, all rights reserved
% [history]: 03/12/12   0.1  new
%-----------------------------------------------------------------------------*/
#include "mex.h"

extern void trop_saast(const double *t, const double *azel, const double *gpos,
					   const double *met_prm, double *tropd, double *tropw);

/* mex interface -------------------------------------------------------------*/
extern void mexFunction(int nargout, mxArray *argout[], int nargin,
					    const mxArray *argin[])
{
	double *t,*azel,*gpos,*met_prm=NULL,*tropd,*tropw=NULL;

    if (nargin<3||mxGetN(argin[0])<1||mxGetN(argin[1])<2||mxGetN(argin[2])<3)
        mexErrMsgTxt("argin error");
	if (nargout>2) mexErrMsgTxt("argout error"); 
    
    t   =mxGetPr(argin[0]);
    azel=mxGetPr(argin[1]);
    gpos=mxGetPr(argin[2]);
    if (nargin>=4&&mxGetN(argin[3])>=3) met_prm=mxGetPr(argin[3]);
	argout[0]=mxCreateDoubleMatrix(1,1,mxREAL); tropd=mxGetPr(argout[0]);
    if (nargout>=2) {
		argout[1]=mxCreateDoubleMatrix(1,1,mxREAL); tropw=mxGetPr(argout[1]);
	}
	trop_saast(t,azel,gpos,met_prm,tropd,tropw);
}
