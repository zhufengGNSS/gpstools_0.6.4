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
% [argout] : posr = station postion (m) (ecef)
%            cdtr = receiver clock bias (m)
% [note]   :
% [version]: $Revision: 15 $ $Date: 06/07/08 1:01 $
%            Copyright(c) 2004-2006 by T.Takasu, all rights reserved
% [history]: 04/10/19  0.1  new
%            05/06/23  0.2  add station index in iz
%            05/08/14  0.3  add argin dflg
%            05/12/15  0.4  set argin posr,cdtr to optionals
%-----------------------------------------------------------------------------*/
#include "mex.h"
#include "qtcmn.h"

extern void PointPos(const double *td, const double *ts, const double *z,
					 const double *iz, int nz, const double *nav,
					 const double *inav, int nnav, const double *posr0,
					 const double *cdtr0, const double *dflg, double *posr,
					 double *cdtr);

/* mex interface -------------------------------------------------------------*/
extern void mexFunction(int nargout, mxArray *argout[], int nargin,
					    const mxArray *argin[])
{
	int nz,nnav;
	double zero[]={0.0,0.0,0.0},df=0.0,c[3];
	double *td,*ts,*z,*iz,*nav,*inav,*posr0,*cdtr0,*dflg,*posr,*cdtr=c;

	if (nargin<6||mxGetM(argin[0])<1||mxGetM(argin[1])<1||
		mxGetM(argin[2])!=mxGetM(argin[3])||mxGetN(argin[3])<2||mxGetN(argin[4])<37||
		mxGetM(argin[4])!=mxGetM(argin[5]))
		mexErrMsgTxt("argin error");
	td  =mxGetPr(argin[0]);
	ts  =mxGetPr(argin[1]);
	z   =mxGetPr(argin[2]); nz=mxGetM(argin[2]);
	iz  =mxGetPr(argin[3]);
	nav =mxGetPr(argin[4]); nnav=mxGetM(argin[4]);
	inav=mxGetPr(argin[5]);
	if (nargin>=7) posr0=mxGetPr(argin[6]); else posr0=zero;
	if (nargin>=8) cdtr0=mxGetPr(argin[7]); else cdtr0=zero;
	if (nargin>=9) dflg=mxGetPr(argin[8]); else dflg=&df;
	
	argout[0]=mxCreateDoubleMatrix(3,1,mxREAL); posr=mxGetPr(argout[0]);
	if (nargout>=2) {
		argout[1]=mxCreateDoubleMatrix(1,1,mxREAL); cdtr=mxGetPr(argout[1]);
	}
	PointPos(td,ts,z,iz,nz,nav,inav,nnav,posr0,cdtr0,dflg,posr,cdtr);
}
