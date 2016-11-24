/*------------------------------------------------------------------------------
% [system] : GpsTools
% [module] : tropospheric mapping function - NMF
% [func]   : calculate tropospheric mapping function by Niell mapping function
% [argin]  : t    = date/time (mjd)
%            azel = azimath/elevation angle(rad) [az,el]
%            gpos = latitude/longitude/height(deg,m) [lat,lon,h]
%            (metprm) = meteological parameters(no use)
% [argout] : mapfd = dry mapping function
%            mapfw = wet mapping function
% [note]   : reference :
%            A.E.Niell, Global mapping functions for the atmosphere delay at
%            radio wavelengths, 1996
% [version]: $Revision: 16 $ $Date: 2008-12-12 15:49:30 +0900 (金, 12 12 2008) $
%            Copyright(c) 2004-2006 by T.Takasu, all rights reserved
% [history]: 04/05/31   0.1  new
%            05/07/01   0.2  fix interpolation bug in case of lat>75
%            08/12/06   0.3  fix bug on hydrostatic seasonal variation
%-----------------------------------------------------------------------------*/
#include "qtcmn.h"

extern void MjdToCal(const double *mjd, const double *sec, double *dt);
extern void CalToMjd(const double *dt, double *mjd, double *sec, int len_dt);

/* linear interpolation ------------------------------------------------------*/
static double Interp1(const double *x, const double *y, int nx, double xi)
{
	double dx;
	int i;
	if (xi<=x[0]) return y[0]; if (x[nx-1]<=xi) return y[nx-1];
	for (i=0;i<nx-1;i++) if (xi<x[i+1]) break;
	dx=(xi-x[i])/(x[i+1]-x[i]);
	return y[i]*(1-dx)+y[i+1]*dx;
}
/* mapping function ----------------------------------------------------------*/
static double tropmf(double el, double a, double b, double c)
{
	return (1.0+a/(1+b/(1.0+c)))/(sin(el)+(a/(sin(el)+b/(sin(el)+c))));
}
/*  tropospheric mapping function - NMF --------------------------------------*/
extern void mapf_nmf(double *t, const double *azel, const double *gpos,
					 double *mapfd, double *mapfw)
{
	static const double lats[]={15.0,30.0,45.0,60.0,75.0};
	static const double nmfh_ave1[]={
		1.2769934e-3, 1.2683230e-3, 1.2465397e-3, 1.2196049e-3, 1.2045996e-3};
	static const double nmfh_ave2[]={
		2.9153695e-3, 2.9152299e-3, 2.9288445e-3, 2.9022565e-3, 2.9024912e-3};
	static const double nmfh_ave3[]={
		62.610505e-3, 62.837393e-3, 63.721774e-3, 63.824265e-3, 64.258455e-3};
	static const double nmfh_amp1[]={
		0.0000000e-0, 1.2709626e-5, 2.6523662e-5, 3.4000452e-5, 4.1202191e-5};
	static const double nmfh_amp2[]={
		0.0000000e-0, 2.1414979e-5, 3.0160779e-5, 7.2562722e-5, 11.723375e-5};
	static const double nmfh_amp3[]={
		0.0000000e-0, 9.0128400e-5, 4.3497037e-5, 84.795348e-5, 170.37206e-5};
	static const double nmfw1[]={
		5.8021897e-4, 5.6794847e-4, 5.8118019e-4, 5.9727542e-4, 6.1641693e-4};
	static const double nmfw2[]={
		1.4275268e-3, 1.5138625e-3, 1.4572752e-3, 1.5007428e-3, 1.7599082e-3};
	static const double nmfw3[]={
		4.3472961e-2, 4.6729510e-2, 4.3908931e-2, 4.4626982e-2, 5.4736038e-2};
	double aht=2.53e-5, bht=5.49e-3, cht=1.14e-3;
	double dt[6],c,ah,bh,ch,aw,bw,cw,dmdh,mjd,lat;
	
	MjdToCal(t,NULL,dt); dt[1]=1.0; dt[2]=28.0;
	CalToMjd(dt,&mjd,NULL,3);
	c=cos(2.0*M_PI*(t[0]-mjd)/365.25)*(gpos[0]<0.0?-1.0:1.0);
	lat=ABS(gpos[0]);
	
	ah=Interp1(lats,nmfh_ave1,5,lat)-Interp1(lats,nmfh_amp1,5,lat)*c;
	bh=Interp1(lats,nmfh_ave2,5,lat)-Interp1(lats,nmfh_amp2,5,lat)*c;
	ch=Interp1(lats,nmfh_ave3,5,lat)-Interp1(lats,nmfh_amp3,5,lat)*c;
	aw=Interp1(lats,nmfw1,5,lat);
	bw=Interp1(lats,nmfw2,5,lat);
	cw=Interp1(lats,nmfw3,5,lat);
	
	/* height correction term */
	dmdh=1/sin(azel[1])-tropmf(azel[1],aht,bht,cht);
	
	/* mapping functions */
	*mapfd=tropmf(azel[1],ah,bh,ch)+dmdh*gpos[2]/1E3;
	*mapfw=tropmf(azel[1],aw,bw,cw);
}
