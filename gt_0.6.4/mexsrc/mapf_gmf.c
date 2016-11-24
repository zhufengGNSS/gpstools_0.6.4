/*------------------------------------------------------------------------------
% [system] : GpsTools
% [module] : tropospheric mapping function - GMF
% [func]   : calculate tropospheric mapping function by Global mapping function
% [argin]  : t    = date/time (mjd)
%            azel = azimath/elevation angle(rad) [az,el]
%            gpos = latitude/longitude/height(deg,m) [lat,lon,h]
%            (metprm) = meteological parameters(no use)
% [argout] : mapfd = dry mapping function
%            mapfw = wet mapping function
% [note]   : reference :
%            J. Boehm et al., Global Mapping Function (GMF): A new empirical
%            mapping function based on numerical weather model data, Geophys.
%            Res. Lett., Vol. 33, L07304, 2006
%            gmf.f (fortran source code) :
%            http://mars.hg.tuwien.ac.at/~ecmwf1/
% [version]: $Revision: 12 $ $Date: 2008-11-25 10:02:15 +0900 (火, 25 11 2008) $
%            Copyright(c) 2004-2006 by T.Takasu, all rights reserved
% [history]: 06/05/05   0.1  new
%-----------------------------------------------------------------------------*/
#include "mex.h"
#include "qtcmn.h"

extern void GMF(double *dmjd, double *dlat, double *dlon, double *dhgt,
				double *zd, double *gmfh, double *gmfw);

/* mex interface -------------------------------------------------------------*/
extern void mexFunction(int nargout, mxArray *argout[], int nargin,
					    const mxArray *argin[])
{
	double *t,*azel,*gpos,*mapfd,*mapfw,lat,lon,hgt,zd;
	
	if (nargin<3||mxGetM(argin[0])<1||mxGetN(argin[1])<2||mxGetN(argin[2])<3)
		mexErrMsgTxt("argin error");
	
	t   =mxGetPr(argin[0]);
	azel=mxGetPr(argin[1]);
	gpos=mxGetPr(argin[2]);
	argout[0]=mxCreateDoubleMatrix(1,1,mxREAL); mapfd=mxGetPr(argout[0]);
	argout[1]=mxCreateDoubleMatrix(1,1,mxREAL); mapfw=mxGetPr(argout[1]);
	
	lat=gpos[0]*DEG2RAD;
	lon=gpos[1]*DEG2RAD;
	hgt=gpos[2];
	zd=M_PI/2.0-azel[1];
	GMF(t,&lat,&lon,&hgt,&zd,mapfd,mapfw);
}
