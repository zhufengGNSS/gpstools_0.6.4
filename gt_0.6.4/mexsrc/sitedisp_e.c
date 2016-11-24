/*------------------------------------------------------------------------------
% [system] : GpsTools
% [module] : station displacement 
% [func]   : calculate station displacement by earth tides
% [argin]  : tut     = date/time(mjd-utc)
%            posr    = station position (m) (ecef)
%            possun  = sun position (m) (ecef)
%            posmoon = moon position (m) (ecef)
%            odisp   = ocean tide amplitudes (m)
%            ophas   = ocean tide phases (deg)
%            gmst    = Greenwich mean sidereal time (rad)
%            xpyp    = pole offset(xp,yp) (rad)
%            opt     = option flags [solid,oload,polar,perm] (1:on)
% [argout] : dpos    = station displacement (m) (ecef)
% [note]   :
% [version]: $Revision: 12 $ $Date: 2008-11-25 10:02:15 +0900 (火, 25 11 2008) $
%            Copyright(c) 2004-2006 by T.Takasu, all rights reserved
% [history]: 04/10/19  0.1  new
%-----------------------------------------------------------------------------*/
#include "qtcmn.h"
#include "mex.h"

/* solar/lunar tides ---------------------------------------------------------*/
static void PlTide(const double *er, const double *posp, double GMp,
				   double lat, double lon, double *dpos)
{
	const double Re=6378137.0,H3=0.292,L3=0.015;
	double rp,ep[3],latp,lonp,p,K2,K3,a,H2,L2,dp,dr,coslp,cosl;
	COPY(posp,3,1,ep); rp=NORM(ep); ep[0]/=rp; ep[1]/=rp; ep[2]/=rp;
	latp=asin(ep[2]);
	lonp=atan2(ep[1],ep[0]);
	K2=GMp/GME*Re*Re*Re*Re/(rp*rp*rp);
	K3=K2*Re/rp;
	
	/* step1 in phase (degree 2) */
	p=(3.0*sin(lat)*sin(lat)-1.0)/2.0;
	H2=0.6078-0.0006*p;
	L2=0.0847+0.0002*p;
	a=DOT(ep,er);
	dp=K2*3.0*L2*a;
	dr=K2*(H2*(1.5*a*a-0.5)-3.0*L2*a*a);
	
	/* step1 in phase (degree 3) */
	dp+=K3*L3*(7.5*a*a-1.5);
	dr+=K3*(H3*(2.5*a*a*a-1.5*a)-L3*(7.5*a*a-1.5)*a);
	
	/* step1 out-of-phase (only radial) */
	coslp=cos(latp); cosl=cos(lat);
	dr=dr+3.0/4.0*0.0025*K2*sin(2.0*latp)*sin(2.0*lat)*sin(lon-lonp);
	dr=dr+3.0/4.0*0.0022*K2*coslp*coslp*cosl*cosl*sin(2.0*(lon-lonp));

	dpos[0]=dp*ep[0]+dr*er[0];
	dpos[1]=dp*ep[1]+dr*er[1];
	dpos[2]=dp*ep[2]+dr*er[2];
}
/* displacement by solid earth tides -----------------------------------------*/
static void SolidTide(const double *possun, const double *posmoon, double lat,
					  double lon, const double *E, double gmst, int opt,
					  double *dpos)
{
	const double GMs=1.327124E+20,GMm=4.902801E+12;
	double dpos1[3],dpos2[3],dr,dn;
	double sinl=sin(lat), sin2l=sin(2.0*lat);
	double er[]={E[2],E[5],E[8]};
	
	/* step1 : time domain */
	PlTide(er,possun,GMs,lat,lon,dpos1);
	PlTide(er,posmoon,GMm,lat,lon,dpos2);
	
	/* step2 : frequency domain, only K1 radial */
	dr=-0.012*sin2l*sin(gmst+lon);

	dpos[0]=dpos1[0]+dpos2[0]+dr*E[2];
	dpos[1]=dpos1[1]+dpos2[1]+dr*E[5];
	dpos[2]=dpos1[2]+dpos2[2]+dr*E[8];

	/* eliminate permanent deformation */
	if (opt) {
		dr=0.1196*(1.5*sinl*sinl-0.5);
		dn=0.0247*sin2l;
		dpos[0]+=dr*E[2]+dn*E[1];
		dpos[1]+=dr*E[5]+dn*E[4];
		dpos[2]+=dr*E[8]+dn*E[7];
	}
}
/* astronomical arguments ----------------------------------------------------*/
static void AstArgs(double tut, double *angs)
{
	static const double args[11][5]={
	{1.40519E-4, 2.0,-2.0, 0.0, 0.00},	/* M2 */
	{1.45444E-4, 0.0, 0.0, 0.0, 0.00},	/* S2 */
	{1.37880E-4, 2.0,-3.0, 1.0, 0.00},	/* N2 */
	{1.45842E-4, 2.0, 0.0, 0.0, 0.00},	/* K2 */
	{0.72921E-4, 1.0, 0.0, 0.0, 0.25},	/* K1 */
	{0.67598E-4, 1.0,-2.0, 0.0,-0.25},  /* O1 */
	{0.72523E-4,-1.0, 0.0, 0.0,-0.25},	/* P1 */
	{0.64959E-4, 1.0,-3.0, 1.0,-0.25},	/* Q1 */
	{0.53234E-5, 0.0, 2.0, 0.0, 0.00},	/* Mf */
	{0.26392E-5, 0.0, 1.0,-1.0, 0.00},	/* Mm */
	{0.03982E-5, 2.0, 0.0, 0.0, 0.00}};	/* Ssa */
	double td,tt,a[5];
	int i,j;
	td=floor(tut);
	tt=(27392.500528+1.000000035*(td-42412.0))/36525.0;
	a[0]=(tut-td)*86400;
	a[1]=(279.69668+(36000.768930485+3.03E-4*tt)*tt)*DEG2RAD; /* H0 */
	a[2]=(((1.9E-6*tt-0.001133)*tt+481267.88314137)*tt+270.434358)*DEG2RAD; /* S0 */
	a[3]=(((-1.2E-5*tt-0.010325)*tt+4069.0340329577)*tt+334.329653)*DEG2RAD; /* P0 */
	a[4]=2.0*M_PI;
	for (i=0;i<11;i++) {
		for (angs[i]=0.0,j=0;j<5;j++) angs[i]+=args[i][j]*a[j];
		angs[i]=angs[i]-floor(angs[i]/(2.0*M_PI))*2.0*M_PI;
	}
}
/* displacement by ocean loading ---------------------------------------------*/
static void OceanLoad(double tut, const double *E, const double *odisp,
					  const double *ophas, double *dpos)
{
	double angs[11],dp[3],dpp[3],Et[9];
	int i,j;
	if (mxIsNaN(odisp[0])) {
		dpos[0]=dpos[1]=dpos[2]=0.0;
		return;
	}
	AstArgs(tut,angs);
	for (i=0;i<3;i++)
	for (j=0,dp[i]=0.0;j<11;j++)
		dp[i]+=odisp[i*11+j]*cos(angs[j]-ophas[i*11+j]*DEG2RAD);
	dpp[0]=-dp[1]; dpp[1]=-dp[2]; dpp[2]=dp[0];
	Tr(E,Et);
	Mv(Et,dpp,dpos);
}
/* displacement by polar tides -----------------------------------------------*/
static void PolarTide(double lat, double lon, const double *E,
					  const double *xpyp, double *dpos)
{
	double Et[9],dp[3],xp=xpyp[0]/SEC2RAD,yp=xpyp[1]/SEC2RAD;
	dp[0]=9E-3*sin(lat)*(xp*sin(lon)+yp*cos(lon));
	dp[1]=-9E-3*cos(2.0*lat)*(xp*cos(lon)-yp*sin(lon));
	dp[2]=-32E-3*sin(2.0*lat)*(xp*cos(lon)-yp*sin(lon));
	Tr(E,Et); Mv(Et,dp,dpos);
}
/* station displacement ------------------------------------------------------*/
extern void SiteDisp(const double *tutc, const double *posr, const double *possun,
					 const double *posmoon, const double *odisp,
					 const double *ophas, const double *gmst, const double *erp,
					 const double *opt, double *dpos)
{
	double lat,lon,sinp,cosp,sinl,cosl,E[9],dp[3];
	lat=asin(posr[2]/NORM(posr));
	lon=atan2(posr[1],posr[0]);
	sinp=sin(lat); cosp=cos(lat);
	sinl=sin(lon); cosl=cos(lon);
	E[0]=-sinl;      E[3]=cosl;       E[6]=0.0; /* XYZ->ENR */
	E[1]=-sinp*cosl; E[4]=-sinp*sinl; E[7]=cosp;
	E[2]=cosp*cosl;  E[5]=cosp*sinl;  E[8]=sinp;
	
	dpos[0]=dpos[1]=dpos[2]=0.0;
	if ((int)opt[0]) {
		SolidTide(possun,posmoon,lat,lon,E,*gmst,(int)opt[3],dp);
		dpos[0]+=dp[0];
		dpos[1]+=dp[1];
		dpos[2]+=dp[2];
	}
	if ((int)opt[1]) {
		OceanLoad(tutc[0]+erp[2]/86400.0,E,odisp,ophas,dp);
		dpos[0]+=dp[0];
		dpos[1]+=dp[1];
		dpos[2]+=dp[2];
	}
	if ((int)opt[2]) {
		PolarTide(lat,lon,E,erp,dp);
		dpos[0]+=dp[0];
		dpos[1]+=dp[1];
		dpos[2]+=dp[2];
	}
}
