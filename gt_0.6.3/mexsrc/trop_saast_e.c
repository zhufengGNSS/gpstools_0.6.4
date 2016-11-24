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
%                       IERS Conventions 2003 ch.9
% [version]: $Revision: 7 $ $Date: 06/07/08 1:01 $
%            Copyright(c) 2004-2006 by T.Takasu, all rights reserved
% [history]: 03/12/12   0.1  new
%            05/06/15   0.2  use iers conventions 2003 for hydrostatic delay
%            05/08/08   0.3  return 0 if height<-1km||32km<height
%-----------------------------------------------------------------------------*/
#include "mex.h"
#include "qtcmn.h"

/* correction terms ----------------------------------------------------------*/
static const double hb[]={0,500,1000,1500,2000,2500,3000,4000,5000};
static const double Bc[]={1.156,1.079,1.006,0.938,0.874,0.813,0.757,0.654,0.563};
static const double hr[]={0,500,1000,1500,2000,3000,4000,5000};
static const double zr[]={0,60,66,70,73,75,76,77,78,78.5,79,79.5,79.75,80,90};
static const double Rc[]={
	0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,
	0.003,0.003,0.002,0.002,0.002,0.002,0.001,0.001,
	0.006,0.006,0.005,0.005,0.004,0.003,0.003,0.002,
	0.012,0.011,0.010,0.009,0.008,0.006,0.005,0.004,
	0.020,0.018,0.017,0.015,0.013,0.011,0.009,0.007,
	0.031,0.028,0.025,0.023,0.021,0.017,0.014,0.011,
	0.039,0.035,0.032,0.029,0.026,0.021,0.017,0.014,
	0.050,0.045,0.041,0.037,0.033,0.027,0.022,0.018,
	0.065,0.059,0.054,0.049,0.044,0.036,0.030,0.024,
	0.075,0.068,0.062,0.056,0.051,0.042,0.034,0.028,
	0.087,0.079,0.072,0.065,0.059,0.049,0.040,0.033,
	0.102,0.093,0.085,0.077,0.070,0.058,0.047,0.039,
	0.111,0.101,0.092,0.083,0.076,0.063,0.052,0.043,
	0.121,0.110,0.100,0.091,0.083,0.068,0.056,0.047,
	0.121,0.110,0.100,0.091,0.083,0.068,0.056,0.047
};
/* linear interpolation ------------------------------------------------------*/
static double Interp1(const double *x, const double *y, int nx, int ny,
					  double xi)
{
	int i;
	double dx;
    
    for (i=0;i<nx-1;i++) if (xi<x[i+1]) break;
    dx=(xi-x[i])/(x[i+1]-x[i]);
    
    return y[i]*(1.0-dx)+y[i+1]*dx;
}
static double Interp2(const double *x, const double *y, const double *z, int nx,
					  int ny, double xi, double yi)
{
	int i,j;
	double dx,dy;
    
    for (i=0;i<nx-1;i++) if (xi<x[i+1]) break;
    for (j=0;j<ny-1;j++) if (yi<y[j+1]) break;
    dx=(xi-x[i])/(x[i+1]-x[i]); dy=(yi-y[j])/(y[j+1]-y[j]);
    
    return z[i*ny+j]*(1.0-dx)*(1.0-dy)+z[(i+1)*ny+j]*dx*(1.0-dy)+
           z[i*ny+j+1]*(1.0-dx)*dy+z[(i+1)*ny+j+1]*dx*dy;
}

/* tropospheric model - Saastamoinen -----------------------------------------*/
extern void trop_saast(const double *t, const double *azel, const double *gpos,
					   const double *met_prm, double *tropd, double *tropw)
{
	const double *met_prm_v;
	double T,z,e,B,dR,tanz,met_prm_d[3],h;
	
	h=gpos[2]<0.0?0.0:gpos[2];
	
	/* return 0 if H>32km */
	if (32000.0<h) {
		*tropd=0.0;
		if (tropw!=NULL) *tropw=0.0;
		return;
	}
	/* standard atmosphere */
    if (met_prm==NULL) {
		met_prm_d[0]=1013.25*pow(1-2.2557E-5*h,5.2568);
		met_prm_d[1]=15.0-6.5E-3*h;
		met_prm_d[2]=50.0;
		met_prm_v=(const double *)met_prm_d;
	}
	else met_prm_v=met_prm;

	/* temperture(K)/zenith angle(rad) */
	T=met_prm_v[1]+273.16;
	z=M_PI/2.0-azel[1]; tanz=tan(z);
	
	/* partial pressure of water vapor(hPa) */
	e=6.108*(met_prm_v[2]/100.0)*exp((17.15*T-4684.0)/(T-38.45));
	
	/* dry tropospheric delay */
	*tropd=0.0022768*met_prm_v[0]/
		   (1.0-0.00266*cos(2.0*gpos[0]*DEG2RAD)-0.00028*h/1E3)/cos(z);
	
	/* wet tropospheric delay */
	if (tropw!=NULL) *tropw=0.002277/cos(z)*(1255.0/T+0.05)*e;
	else *tropd+=0.002277/cos(z)*(1255.0/T+0.05)*e;
}
