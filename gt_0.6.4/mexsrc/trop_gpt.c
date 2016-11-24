/*------------------------------------------------------------------------------
% [system] : GpsTools
% [module] : tropospheric model - GPT
% [func]   : global pressure and temperature based spherical harmonics
% [argin]  : t    = date/time (mjd)
%            gpos = latitude/longitude/height(deg,m) [lat,lon,h]
% [argout] : pres = pressure (hPa)
%            temp = temperature (C)
%            undu = geoid undulation (m)
% [note]   : reference :
%            J. Boehm, R. Heinkelmann, and H. Schuh (2007), Short Note: A global
%            model of pressure and temperature for geodetic applications,
%            Journal of Geodesy , doi:10.1007/s00190-007-0135-3
%            gpt.f (fortran source code) :
%            http://mars.hg.tuwien.ac.at/~ecmwf1/
% [version]: $Revision: 16 $ $Date: 2008-12-12 15:49:30 +0900 (金, 12 12 2008) $
%            Copyright(c) 2004-2006 by T.Takasu, all rights reserved
% [history]: 08/11/24   0.1  new
%-----------------------------------------------------------------------------*/
#include "mex.h"
#include "qtcmn.h"

extern void GPT(double *dmjd, double *dlat, double *dlon, double *dhgt,
				double *pres, double *temp, double *undu);

/* mex interface -------------------------------------------------------------*/
extern void mexFunction(int nargout, mxArray *argout[], int nargin,
					    const mxArray *argin[])
{
	double *t,*gpos,lat,lon,hgt,pres,temp,undu;
	
	if (nargin<2||mxGetM(argin[0])<1||mxGetN(argin[1])<3)
		mexErrMsgTxt("argin error");
	
	t=mxGetPr(argin[0]);
	gpos=mxGetPr(argin[1]);
	lat=gpos[0]*DEG2RAD;
	lon=gpos[1]*DEG2RAD;
	hgt=gpos[2];
	
	GPT(t,&lat,&lon,&hgt,&pres,&temp,&undu);
	
	if (nargout>0) {
		argout[0]=mxCreateDoubleMatrix(1,1,mxREAL); *mxGetPr(argout[0])=pres;
	}
	if (nargout>1) {
		argout[1]=mxCreateDoubleMatrix(1,1,mxREAL); *mxGetPr(argout[1])=temp;
	}
	if (nargout>2) {
		argout[2]=mxCreateDoubleMatrix(1,1,mxREAL); *mxGetPr(argout[2])=undu;
	}
}
