/*------------------------------------------------------------------------------
% [system] : GpsTools
% [module] : ionospheric delay model - Klobuchar model
% [func]   : calculate ionospheric delay(L1) by Kloubachar model
% [argin]  : t    = date/time(mjd-utc)
%            azel = satellite azimath/elevation angle(rad) [az,el]
%            gpos = latitude/longitude/height(deg,m) [lat,lon,h]
%            (ion_prm) = ionospheric parameters [alpha0,1,2,3;beta0,1,2,3]
%            <global>
%            utc_tai = utc-tai(sec) (default:-32)
% [argout] : ion  = ionospheric delay(L1)(m)
% [note]   : Reference: ICD-GPS-200C Navstar GPS Space Segment/Navigation User
%                       Interfaces Rev.C, 20.3.3.5.2.5
% [version]: $Revision: 16 $ $Date: 2008-12-12 15:49:30 +0900 (金, 12 12 2008) $
%            Copyright(c) 2004-2006 by T.Takasu, all rights reserved
% [history]: 03/12/31   0.1  new
%            08/11/20   0.2  mexGetArray -> mexGetVariable (gt_0.6.4)
%-----------------------------------------------------------------------------*/
#include "mex.h"
#include "qtcmn.h"

static double *utc_tai=NULL;

/* ionospheric model - Klobuchar model ---------------------------------------*/
extern double ion_klob(const double *t, const double *azel, const double *gpos,
					   const double *ion_prm)
{
	static const double ion_prm_d[]={ /* default parameters */
		/* alpha       beta  */
		 0.1676E-7,  0.1352E+6,
		-0.7451E-8, -0.1802E+6,
		-0.5960E-7,  0.6554E+5,
		 0.1192E-6,  0.0000E+0
	};
	static double utc_tai_d[]={-32.0};
	double F,psi,phii,lami,phim,tt,AMP,PER,x;
	
	if (ion_prm==NULL) ion_prm=ion_prm_d;
	if (utc_tai==NULL) utc_tai=utc_tai_d;

	/* earth centered angle(semi-circle) */
	psi=0.0137/(azel[1]/M_PI+0.11)-0.022;
	
	/* subionospheric latitude/longitude(semi-circle) */
	phii=gpos[0]/180.0+psi*cos(azel[0]);
	if (phii<-0.416) phii=-0.416; else if (phii>0.416) phii=0.416;
	lami=gpos[1]/180.0+psi*sin(azel[0])/cos(phii*M_PI);
	
	/* geomagnetic latitude(semi-circle) */
	phim=phii+0.064*cos((lami-1.617)*M_PI);
	
	/* local time(sec) */
	tt=4.32E4*lami+(t[0]-floor(t[0]))*86400.0-utc_tai[0]-19.0;
	if (tt<0) tt=tt+86400.0; else if (tt>=86400.0) tt=tt-86400.0;
	
	/* slant factor */
	F=1.0+16.0*pow(0.53-azel[1]/M_PI,3.0);
	
	/* ionospheric delay */
	AMP=ion_prm[0]+phim*(ion_prm[2]+phim*(ion_prm[4]+phim*ion_prm[6]));
	PER=ion_prm[1]+phim*(ion_prm[3]+phim*(ion_prm[5]+phim*ion_prm[7]));
	AMP=MAX(AMP,0.0);
	PER=MAX(PER,72000.0);
	
	x=2.0*M_PI*(tt-50400.0)/PER;
	if (ABS(x)<1.57) return M_C*F*(5E-9+AMP*(1.0+x*x*(-0.5+x*x/24.0)));
	return M_C*F*5E-9;
}
/* mex interface -------------------------------------------------------------*/
extern void mexFunction(int nargout, mxArray *argout[], int nargin,
					    const mxArray *argin[])
{
	mxArray *p;
	double *t,*azel,*gpos,*ion_prm=NULL;

    if (nargin<3||mxGetN(argin[0])<1||mxGetN(argin[1])<2||mxGetN(argin[2])<3)
        mexErrMsgTxt("argin error");
	if (nargout>1) mexErrMsgTxt("argout error"); 
    
    t=mxGetPr(argin[0]);
    azel=mxGetPr(argin[1]);
    gpos=mxGetPr(argin[2]);
    if (nargin>=4&&mxGetN(argin[3])>=4&&mxGetM(argin[3])>=2) ion_prm=mxGetPr(argin[3]);
	if ((p=mexGetVariable("global","utc_tai"))&&mxGetM(p)>=1) utc_tai=mxGetPr(p);
	argout[0]=mxCreateDoubleMatrix(1,1,mxREAL);
    *mxGetPr(argout[0])=ion_klob(t,azel,gpos,ion_prm);
}
