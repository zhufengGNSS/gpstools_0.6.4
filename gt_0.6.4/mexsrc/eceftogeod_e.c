/*------------------------------------------------------------------------------
% [system] : GpsTools
% [module] : ecef to geodetic position
% [func]   : convert ecef postion to geodetic (latitude,longitude,height)
% [argin]  : epos = ecef position (m) [x;y;z]
% [argout] : gpos = latitude/longitude/height(deg,m) [lat,lon,h] 
%           (E)   = ecef to local tangental coordinate transformation matrix
% [note]   : WGS84 ellipsoide
% [version]: $Revision: 12 $ $Date: 2008-11-25 10:02:15 +0900 (火, 25 11 2008) $
%            Copyright(c) 2004-2006 by T.Takasu, all rights reserved
% [history]: 04/05/31   0.1  new
%            05/12/12   0.2  add argout E
%-----------------------------------------------------------------------------*/
#include "qtcmn.h"

/* ecef to geodetic position -------------------------------------------------*/
extern void EcefToGeod(const double *epos, double *gpos, double *E)
{
	static const double Re=6378137.0;			/* WGS84 */
	static const double f=1.0/298.257223563;	/* WGS84 */
	double e2,r2,z,zz,n,sinp,cosp,sinl,cosl,lat;

	e2=f*(2.0-f);
	r2=epos[0]*epos[0]+epos[1]*epos[1];
	z=epos[2]; zz=0; n=Re;
	while (ABS(z-zz)>=1E-4) {
	    zz=z;
	    sinp=z/sqrt(r2+z*z);
	    n=Re/sqrt(1.0-e2*sinp*sinp);
	    z=epos[2]+n*e2*sinp;
	}
	if (r2<1E-12) lat=90.0; else lat=atan(z/sqrt(r2))*180.0/M_PI;
	gpos[0]=lat;
	gpos[1]=atan2(epos[1],epos[0])*180.0/M_PI;
	gpos[2]=sqrt(r2+z*z)-n;
	
	/* ecef to local tangental coordinate transformation matrix */
	if (E!=NULL) {
		sinp=sin(gpos[0]*DEG2RAD); cosp=cos(gpos[0]*DEG2RAD);
		sinl=sin(gpos[1]*DEG2RAD); cosl=cos(gpos[1]*DEG2RAD);
		E[0]=-sinl;      E[3]=cosl;       E[6]=0.0;
		E[1]=-sinp*cosl; E[4]=-sinp*sinl; E[7]=cosp;
		E[2]=cosp*cosl;  E[5]=cosp*sinl;  E[8]=sinp;
	}
}
