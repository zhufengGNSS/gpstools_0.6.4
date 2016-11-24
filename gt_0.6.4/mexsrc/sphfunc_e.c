/*------------------------------------------------------------------------------
% [system] : GpsTools
% [module] : spheric harmonic functions
% [func]   : spheric harmonic functions by azimuth/elevation angle
% [argin]  : az,el = azimuth/elevation angle (rad)
%            nmax  = max degrees of spheric harmonic functions
% [argout] : fc,fs = spheric harmonic functions
%                fc(n+1,m+1) = Pnm(-cos(2*el))*cos(m*az)
%                fs(n+1,m+1) = Pnm(-cos(2*el))*sin(m*az)
%                (Pnm = normalized Legendre polynomial)
% [note]   :
% [version]: $Revision: 12 $ $Date: 2008-11-25 10:02:15 +0900 (火, 25 11 2008) $
%            Copyright(c) 2004-2006 by T.Takasu, all rights reserved
% [history]: 04/10/18  0.1  new
%-----------------------------------------------------------------------------*/
#include "qtcmn.h"

#define NM(n,m)		((n)+(m)*(nmax+1))		/* pnm,dpnm index */

/* factorial -----------------------------------------------------------------*/
static double factorial(int n)
{
	double f; for (f=1.0;n>1;n--) f*=(double)n; return f;
}
/* normalized legendre functions ---------------------------------------------*/
static void LegendreFunc(double x, int nmax, double *p)
{
	int n,m;
	
	p[NM(0,0)]=1.0;
	p[NM(1,0)]=x;
	p[NM(1,1)]=sqrt(1.0-x*x);
	
	for (n=2;n<=nmax;n++) {
	    p[NM(n,n)]=(2.0*n-1.0)*sqrt(1.0-x*x)*p[NM(n-1,n-1)];
	    for (m=0;m<n;m++) {
	        p[NM(n,m)]=(x*(2.0*n-1.0)*p[NM(n-1,m)]-(n+m-1.0)*p[NM(n-2,m)])/(n-m);
	    }
	}
	for (n=1;n<=nmax;n++) {
	    p[NM(n,0)]*=sqrt(2.0*n+1.0);
	    for (m=1;m<=n;m++)
	    	p[NM(n,m)]*=sqrt(factorial(n-m)*(4.0*n+2.0)/factorial(n+m));
	}
}
/* spheric harmonic functions ------------------------------------------------*/
extern void SphFunc(double az, double el, int nmax, double *fc, double *fs)
{
	double *p,*cosm,*sinm;
	int n,m;
	
	p=ZMAT(nmax+1,nmax+1); cosm=VEC(nmax+1); sinm=VEC(nmax+1);
	if (p==NULL||cosm==NULL||sinm==NULL) return;
	
	for (m=0;m<=nmax;m++) {
		cosm[m]=cos((double)m*az);
		sinm[m]=sin((double)m*az);
	}
	LegendreFunc(-cos(2.0*el),nmax,p);
	
	fc[0]=1.0;
	for (n=1;n<=nmax;n++)
	for (m=0;m<=n;m++) {
		fc[NM(n,m)]=p[NM(n,m)]*cosm[m];
		fs[NM(n,m)]=p[NM(n,m)]*sinm[m];
	}
	FreeMat(p); FreeMat(cosm); FreeMat(sinm);
}
