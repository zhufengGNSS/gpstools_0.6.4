/*------------------------------------------------------------------------------
% [system] : GpsTools
% [module] : eci to satellite-fixed coordinate transformation matrix
% [func]   : calculate eci to satellite-fixed coordinate transformation matrix
% [argin]  : state = satellite position/velocity[x;y;z;vx;vy;vz](m) (eci)
% [argout] : E = eci to satellite-fixed coordinate transformation matrix(3x3)
%                [r_radial;r_along-track;r_cross-track]=E*r
% [note]   :
% [version]: $Revision: 12 $ $Date: 2008-11-25 10:02:15 +0900 (火, 25 11 2008) $
%            Copyright(c) 2004-2006 by T.Takasu, all rights reserved
% [history]: 04/02/08   0.1  new
%            04/08/11   0.2  separated from ecsftosatf.c
%-----------------------------------------------------------------------------*/
#include "qtcmn.h"

/* eci to satellite-fixed coordinate transformation matrix -------------------*/
extern void EcsfToSatf(const double *state, double *E)
{
	double crt[3],alt[3],nrad,ncrt,nalt;
	
	CROSS(state,state+3,crt);
	CROSS(crt,state,alt);
	nrad=NORM(state);
	nalt=NORM(alt);
	ncrt=NORM(crt);
	E[0]=state[0]/nrad;
	E[3]=state[1]/nrad;
	E[6]=state[2]/nrad;
	E[1]=alt[0]/nalt;
	E[4]=alt[1]/nalt;
	E[7]=alt[2]/nalt;
	E[2]=crt[0]/ncrt;
	E[5]=crt[1]/ncrt;
	E[8]=crt[2]/ncrt;
}
