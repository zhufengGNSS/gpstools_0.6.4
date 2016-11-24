/*------------------------------------------------------------------------------
% [system] : GpsTools
% [module] : shadow function
% [func]   : calculate shadow factor in eclipse
% [argin]  : rpos  = satellite position (m)
%            rsun  = sun position (m)
%            rmoon = moon position (m)
%           (model)= shadow model
%                    'cylind'  : cylindric model
%                    'penumbra': penumbra/umbra model
% [argout] : fact  = shadow factor (0:in eclipse,1:not in eclipse)
% [note]   :
% [version]: $Revision: 4 $ $Date: 06/07/08 1:01 $
%            Copyright(c) 2004-2006 by T.Takasu, all rights reserved
% [history]: 04/09/28  0.1  new
%            05/04/01  0.2  add penumbra/umbra model
%-----------------------------------------------------------------------------*/
#include <string.h>
#include "qtcmn.h"

/* penumbra/umbra shadow model -----------------------------------------------*/
static double pshadow(const double *rpos, const double *rsun)
{
	static const double Re=6378136.3;
	static const double Rs=6.960E+8;
	int i;
	double rrs[3],f,fe,fs,pe,ps;
	for (i=0;i<3;i++) rrs[i]=rsun[i]-rpos[i];
	f=M_PI-acos(DOT(rpos,rrs)/NORM(rpos)/NORM(rrs));
	fe=asin(Re/NORM(rpos));
	fs=asin(Rs/NORM(rrs));
	if (f>=fe+fs) return 1.0; /* sunlight */
	if (f<=fe-fs) return 0.0; /* umbra */
	
	/* penumbra */
	pe=acos((f*f+fe*fe-fs*fs)/(2*f*fe));
	ps=acos((f*f+fs*fs-fe*fe)/(2*f*fs));
	return 1.0-(fe*fe*pe+fs*fs*ps-f*fe*sin(pe))/(M_PI*fs*fs);
}
/* shadow function -----------------------------------------------------------*/
extern double ShadowFunc(const double *rpos, const double *rsun,
					     const double *rmoon, const char *model)
{
	static const double Re=6378136.3;
	static const double Rm=1738000.0;
	int i;
	double rrm[3],rsm[3],phi,fact=1.0;
	
	if (strcmp(model,"penumbra")==0) {
		fact=pshadow(rpos,rsun);
	}
	else {
		phi=acos(DOT(rpos,rsun)/NORM(rpos)/NORM(rsun));
		if (phi>=M_PI/2.0&&NORM(rpos)*sin(phi)<=Re) return 0.0; /* earth shadow */
	}
	for (i=0;i<3;i++) {
		rrm[i]=rpos[i]-rmoon[i];
		rsm[i]=rsun[i]-rmoon[i];
	}
	phi=acos(DOT(rrm,rsm)/NORM(rrm)/NORM(rsm));
	if (phi>=M_PI/2.0&&NORM(rrm)*sin(phi)<=Rm) return 0.0; /* moon shadow */
	return fact;
}
