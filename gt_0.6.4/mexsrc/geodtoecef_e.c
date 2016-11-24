/*------------------------------------------------------------------------------
% [system] : GpsTools
% [module] : geodetic to ecef position
% [func]   : convert geodetic postion to ecef
% [argin]  : gpos = latitude/longitude/height(deg,m) [lat,lon,h] 
% [argout] : epos = ecef position (m) [x;y;z]
%            E    = ecef to local tangental coordinate transformation matrix
% [note]   : WGS84 ellipsoide
% [version]: $Revision: 12 $ $Date: 2008-11-25 10:02:15 +0900 (火, 25 11 2008) $
%            Copyright(c) 2004-2006 by T.Takasu, all rights reserved
% [history]: 04/05/31   0.1  new
%-----------------------------------------------------------------------------*/
#include "qtcmn.h"

/* geodetic to ecef position -------------------------------------------------*/
extern void GeodToEcef(const double *gpos, double *epos, double *E)
{
	static const double Re=6378137.0;			/* WGS84 */
	static const double f=1.0/298.257223563;	/* WGS84 */
	double e2,sinp,cosp,sinl,cosl,n;

	e2=f*(2.0-f);
	sinp=sin(gpos[0]*DEG2RAD); cosp=cos(gpos[0]*DEG2RAD);
	sinl=sin(gpos[1]*DEG2RAD); cosl=cos(gpos[1]*DEG2RAD);
	n=Re/sqrt(1.0-e2*sinp*sinp);
	epos[0]=(n+gpos[2])*cosp*cosl;
	epos[1]=(n+gpos[2])*cosp*sinl;
	epos[2]=((1.0-e2)*n+gpos[2])*sinp;

	/* ecef to local tangental coordinate transformation matrix */
	if (E!=NULL) {
	    E[0]=-sinl;      E[3]=cosl;       E[6]=0.0;
	    E[1]=-sinp*cosl; E[4]=-sinp*sinl; E[7]=cosp;
	    E[2]=cosp*cosl;  E[5]=cosp*sinl;  E[8]=sinp;
	}
}
