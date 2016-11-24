/*-------------------------------------------------------------------------------
% [system] : GpsTools
% [module] : satellite azimuth/elevation angle
% [func]   : calculate satellite azimuth/elevation angle
% [argin]  : spos = satellite postion [posx;posy;posz] (m) (ecef)
%            rpos = station position [posx;posy;posz] (m) (ecef)
% [argout] : azel = azimuth/elevation angle [az,el] (rad)
% [note]   :
% [version]: $Revision: 16 $ $Date: 2008-12-12 15:49:30 +0900 (金, 12 12 2008) $
%            Copyright(c) 2004-2008 by T.Takasu, all rights reserved
% [history]: 04/05/31   0.1  new
%            08/11/30   0.2  suppress warning
%-----------------------------------------------------------------------------*/
#include "qtcmn.h"

extern void EcefToGeod(const double *epos, double *gpos, double *E);

/* satellite azimuth/elevation angle -----------------------------------------*/
extern void SatAzEl(const double *spos, const double *rpos, double *azel)
{
	double gpos[3],rs[3],rss[3],E[9];
	EcefToGeod(rpos,gpos,E);
	rs[0]=spos[0]-rpos[0];
	rs[1]=spos[1]-rpos[1];
	rs[2]=spos[2]-rpos[2];
	Mv(E,rs,rss);
	azel[0]=atan2(rss[0],rss[1]);
	azel[1]=asin(rss[2]/NORM(rs));
}
