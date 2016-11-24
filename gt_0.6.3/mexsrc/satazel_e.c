/*-------------------------------------------------------------------------------
% [system] : GpsTools
% [module] : satellite azimuth/elevation angle
% [func]   : calculate satellite azimuth/elevation angle
% [argin]  : spos = satellite postion [posx;posy;posz] (m) (ecef)
%            rpos = station position [posx;posy;posz] (m) (ecef)
% [argout] : azel = azimuth/elevation angle [az,el] (rad)
% [note]   :
% [version]: $Revision: 4 $ $Date: 06/07/08 1:01 $
%            Copyright(c) 2004-2006 by T.Takasu, all rights reserved
% [history]: 04/05/31   0.1  new
%-----------------------------------------------------------------------------*/
#include "qtcmn.h"

extern void GeodToEcef(const double *gpos, double *epos, double *E);

/* satellite azimuth/elevation angle -----------------------------------------*/
extern void SatAzEl(const double *spos, const double *rpos, double *azel)
{
	double gpos[3],epos[3],rs[3],rss[3],E[9];
	EcefToGeod(rpos,gpos,E);
	rs[0]=spos[0]-rpos[0];
	rs[1]=spos[1]-rpos[1];
	rs[2]=spos[2]-rpos[2];
	Mv(E,rs,rss);
	azel[0]=atan2(rss[0],rss[1]);
	azel[1]=asin(rss[2]/NORM(rs));
}
