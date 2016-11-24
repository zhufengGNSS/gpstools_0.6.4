/*------------------------------------------------------------------------------
% [system] : GpsTools
% [module] : tropospheric mapping function - VMF1
% [func]   : calculate tropospheric mapping function by Vienna mapping function 1
% [argin]  : t    = date/time (mjd)
%            azel = azimath/elevation angle(rad) [az,el]
%            gpos = latitude/longitude/height(deg,m) [lat,lon,h]
%            ah   = hydrostatic coefficient a
%            aw   = wet coefficient a
% [argout] : mapfd = dry mapping function
%            mapfw = wet mapping function
% [note]   : reference :
%            J.Boehm et al., Troposphere mapping functions for GPS and very long
%            baseline interferometry from European Centre for Medium-Range
%            Weather Forecasts operational analysis data, J. Geoph. Res.,
%            Vol. 111, B02406, 2006
%            vmf1_ht.f (fortran source code) :
%            http://www.hg.tuwien.ac.at/~ecmwf1/
% [version]: $Revision: 2 $ $Date: 06/07/08 1:01 $
%            Copyright(c) 2004-2008 by T.Takasu, all rights reserved
% [history]: 08/12/06   0.1  new
%-----------------------------------------------------------------------------*/
#include "mex.h"
#include "qtcmn.h"

extern void VMF1_HT(double *ah, double *aw, double *dmjd, double *dlat,
				    double *dhgt, double *zd, double *vmf1h, double *vmf1w);

/* mex interface -------------------------------------------------------------*/
extern void mexFunction(int nargout, mxArray *argout[], int nargin,
					    const mxArray *argin[])
{
	double *t,*azel,*gpos,*ah,*aw,*mapfd,*mapfw,lat,hgt,zd;
	
	if (nargin<5||mxGetM(argin[0])<1||mxGetN(argin[1])<2||mxGetN(argin[2])<3)
		mexErrMsgTxt("argin error");
	
	t   =mxGetPr(argin[0]);
	azel=mxGetPr(argin[1]);
	gpos=mxGetPr(argin[2]);
	ah  =mxGetPr(argin[3]);
	aw  =mxGetPr(argin[4]);
	argout[0]=mxCreateDoubleMatrix(1,1,mxREAL); mapfd=mxGetPr(argout[0]);
	argout[1]=mxCreateDoubleMatrix(1,1,mxREAL); mapfw=mxGetPr(argout[1]);
	
	lat=gpos[0]*DEG2RAD;
	hgt=gpos[2];
	zd=M_PI/2.0-azel[1];
	VMF1_HT(ah,aw,t,&lat,&hgt,&zd,mapfd,mapfw);
}
