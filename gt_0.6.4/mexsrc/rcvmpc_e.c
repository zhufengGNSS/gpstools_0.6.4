/*------------------------------------------------------------------------------
% [system] : GpsTools
% [module] : multipath bias correction
% [func]   : multipath bias correction
% [argin]  : azel = azimuth/elevation angle(rad) [az,el]
%            mpcc = multipath bias cooefficient Cnm
%            mpcs = multipath bias cooefficient Snm
% [argout] : mpr  = multipath bias (m)
% [note]   :
% [version]: $Revision: 12 $ $Date: 2008-11-25 10:02:15 +0900 (火, 25 11 2008) $
%            Copyright(c) 2004-2006 by T.Takasu, all rights reserved
% [history]: 04/10/19   0.1  new
%-----------------------------------------------------------------------------*/
#include "qtcmn.h"
#include "mex.h"

#define NM(n,m)		((n)+(m)*(nmax+1))		/* fs,fc index */

extern void SphFunc(double az, double el, int nmax, double *fc, double *fs);

/* multipath bias correction -------------------------------------------------*/
extern double RcvMpc(const double *azel, const double *mpcc, const double *mpcs,
					 int nmax)
{
	double *fc,*fs,mpr=0.0;
	int n,m;
	if (nmax<=0) return 0.0;
	for (n=0;n<=nmax;n++)
	for (m=0;m<=n;m++) if (mxIsNaN(mpcc[NM(n,m)])) return 0.0;
	
	fs=ZMAT(nmax+1,nmax+1);
	fc=ZMAT(nmax+1,nmax+1);
	if (fs==NULL||fc==NULL) mexErrMsgTxt("malloc error");
	
	SphFunc(azel[0],azel[1],nmax,fc,fs);
	for (n=0;n<=nmax;n++)
	for (m=0;m<=n;m++) {
		mpr+=mpcc[NM(n,m)]*fc[NM(n,m)]+mpcs[NM(n,m)]*fs[NM(n,m)];
	}
	FreeMat(fs); FreeMat(fc);
	
	return mpr;
}
