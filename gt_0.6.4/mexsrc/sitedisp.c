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
#include "mex.h"

extern void SiteDisp(const double *tutc, const double *posr, const double *possun,
					 const double *posmoon, const double *odisp,
					 const double *ophas, const double *gmst, const double *erp,
					 const double *opt, double *dpos);

/* mex interface -------------------------------------------------------------*/
extern void mexFunction(int nargout, mxArray *argout[], int nargin,
					    const mxArray *argin[])
{
	const double opt_v[]={1.0,1.0,1.0,0.0}, *opt;
	double *tutc,*posr,*possun,*posmoon,*odisp,*ophas,*gmst,*erp,*dpos;
	
	if (nargin<8||mxGetM(argin[0])<1||mxGetM(argin[1])<3||mxGetM(argin[2])<3||
		mxGetM(argin[3])<3||mxGetM(argin[4])<11||mxGetN(argin[4])<3||
		mxGetM(argin[5])<11||mxGetN(argin[5])<3||mxGetM(argin[6])<1||
		mxGetN(argin[7])<3)
		mexErrMsgTxt("argin error");
	
	tutc   =mxGetPr(argin[0]);
	posr   =mxGetPr(argin[1]);
	possun =mxGetPr(argin[2]);
	posmoon=mxGetPr(argin[3]);
	odisp  =mxGetPr(argin[4]);
	ophas  =mxGetPr(argin[5]);
	gmst   =mxGetPr(argin[6]);
	erp    =mxGetPr(argin[7]);
	if (nargin>8&&mxGetN(argin[8])>=4) opt=mxGetPr(argin[8]); else opt=opt_v;

	argout[0]=mxCreateDoubleMatrix(3,1,mxREAL); dpos=mxGetPr(argout[0]);
	
	SiteDisp(tutc,posr,possun,posmoon,odisp,ophas,gmst,erp,opt,dpos);
}
