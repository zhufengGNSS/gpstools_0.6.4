/*-------------------------------------------------------------------------------
% [system] : GpsTools
% [module] : point positioning
% [argin]  : td   = date(mjd-gpst)
%            ts   = time(sec)
%            z    = observation data(pseudo-range)
%            iz   = observation data satellite/station index
%            nav  = navigation messages
%            inav = navigation messages satellite index
%           (posr)= approx station position (m) (ecef)
%           (cdtr)= approx receiver clock bias (m)
%           (dflg)= debug flag (0:off,1:on,2:detailed)
%           (opt) = option (0:none,1:tropos delay off)
% [argout] : posr = station postion (m) (ecef)
%            cdtr = receiver clock bias (m)
% [note]   :
% [version]: $Revision: 16 $ $Date: 2008-12-12 15:49:30 +0900 (金, 12 12 2008) $
%            Copyright(c) 2004-2008 by T.Takasu, all rights reserved
% [history]: 04/10/19  0.1  new
%            05/06/23  0.2  add station index in iz
%            05/08/14  0.3  add argin dflg
%            05/12/15  0.4  set argin posr,cdtr to optionals
%            08/11/21  0.5  add argin opt (gt_0.6.4)
%                           suppress warning
%-----------------------------------------------------------------------------*/
#include "mex.h"
#include "qtcmn.h"

extern void PointPos(const double *td, const double *ts, const double *z,
					 const double *iz, int nz, const double *nav,
					 const double *inav, int nnav, const double *posr0,
					 const double *cdtr0, const double *dflg, const double *opt,
					 double *posr, double *cdtr);

/* mex interface -------------------------------------------------------------*/
extern void mexFunction(int nargout, mxArray *argout[], int nargin,
					    const mxArray *argin[])
{
	int nz,nnav;
	double zero[]={0.0,0.0,0.0},df=0.0,c[3];
	double *td,*ts,*z,*iz,*nav,*inav,*posr0,*cdtr0,*dflg,*posr,*opt,*cdtr=c;

	if (nargin<6||mxGetM(argin[0])<1||mxGetM(argin[1])<1||
		mxGetM(argin[2])!=mxGetM(argin[3])||mxGetN(argin[3])<2||mxGetN(argin[4])<37||
		mxGetM(argin[4])!=mxGetM(argin[5]))
		mexErrMsgTxt("argin error");
	td  =mxGetPr(argin[0]);
	ts  =mxGetPr(argin[1]);
	z   =mxGetPr(argin[2]); nz=(int)mxGetM(argin[2]);
	iz  =mxGetPr(argin[3]);
	nav =mxGetPr(argin[4]); nnav=(int)mxGetM(argin[4]);
	inav=mxGetPr(argin[5]);
	if (nargin>=7) posr0=mxGetPr(argin[6]); else posr0=zero;
	if (nargin>=8) cdtr0=mxGetPr(argin[7]); else cdtr0=zero;
	if (nargin>=9) dflg=mxGetPr(argin[8]); else dflg=&df;
	if (nargin>=10) opt=mxGetPr(argin[0]); else opt=zero;
	
	argout[0]=mxCreateDoubleMatrix(3,1,mxREAL); posr=mxGetPr(argout[0]);
	if (nargout>=2) {
		argout[1]=mxCreateDoubleMatrix(1,1,mxREAL); cdtr=mxGetPr(argout[1]);
	}
	PointPos(td,ts,z,iz,nz,nav,inav,nnav,posr0,cdtr0,dflg,opt,posr,cdtr);
}
