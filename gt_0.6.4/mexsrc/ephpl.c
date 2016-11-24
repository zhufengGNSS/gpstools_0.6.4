/*------------------------------------------------------------------------------
% [system] : GpsTools
% [module] : solar/planetary ephemeris
% [func]   : read solar/planetary postions from de405 ephemeris
% [argin]  : t = date/time(mjd-utc)
%            pl= solar/planetary numbers
%                (1:mercury,2:venus,3:earth-moon barycenter,4:mars,5:jupiter,
%                 6:saturn,7:uranus,8:neptune,9:pluto,10:moon,11:sun)
%            <global>
%            utc_tai = utc-tai(sec)
%            ephpdir = ephemris data directory
% [argout] : r = solar/planetary positions(km)
% [note]   : coordinate:barycenter(all except for moon),eci(icrf)(moon)
% [version]: $Revision: 16 $ $Date: 2008-12-12 15:49:30 +0900 (金, 12 12 2008) $
%            Copyright(c) 2004-2008 by T.Takasu, all rights reserved
% [history]: 04/03/13  0.1  new
%            08/11/20  0.2  mexGetArray -> mexGetVariable (gt_0.6.4)
%                           suppress warning
%-----------------------------------------------------------------------------*/
#include "mex.h"

static double utc_tai_v=-32.0;
double *utc_tai=&utc_tai_v;
extern char ephpdir[];
extern int EphPl(const double *t, const double *pl, int np, double *r);

/* mex interface -------------------------------------------------------------*/
extern void mexFunction(int nargout, mxArray *argout[], int nargin,
					    const mxArray *argin[])
{
	mxArray *p;
	int np;
	double *t,*pl,*r;

	if (nargin<2||mxGetM(argin[0])<1||mxGetM(argin[1])<1)
		mexErrMsgTxt("argin error");
	
	t=mxGetPr(argin[0]);
	pl=mxGetPr(argin[1]);
	if (mxGetM(argin[1])==1) np=(int)mxGetN(argin[1]); else np=(int)mxGetM(argin[1]);
	if ((p=mexGetVariable("global","utc_tai"))&&mxGetM(p)>=1) utc_tai=mxGetPr(p);
	if ((p=mexGetVariable("global","ephpdir"))&&mxIsChar(p))
		mxGetString(p,ephpdir,256);
		
	argout[0]=mxCreateDoubleMatrix(3,np,mxREAL); r=mxGetPr(argout[0]);
	if (!EphPl(t,pl,np,r)) mexErrMsgTxt("planet ephemeris read error");
}
