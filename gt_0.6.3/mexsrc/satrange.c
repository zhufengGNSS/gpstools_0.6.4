/*-------------------------------------------------------------------------------
% [system] : GpsTools
% [module] : satellite-station range
% [func]   : calculate range between satellite and station
% [argin]  : state = satellite position/velocity (m) (eci)
%            rpos  = station position (m) (ecef)
%            U     = eci to ecef transformation matrix
%            dtr   = receiver clock bias (sec)
%            (corrlightt) = light-time/receiver clock correction flag (1:on)
%            (direc) = ranging direction (1:sat to sta,2:sta to sat)
% [argout] : rs    = satellite tx position(m) (eci)
%            range = range between satellite and station(m)
%            (drds) = partical derivatives of range by satellite state
% [note]   :
% [version]: $Revision: 4 $ $Date: 06/07/08 1:01 $
%            Copyright(c) 2004-2006 by T.Takasu, all rights reserved
% [history]: 04/04/12   0.1  new
%-----------------------------------------------------------------------------*/
#include "mex.h"

extern void SatRange(const double *state, const double *rpos, const double *U,
					 const double *dtr, int corrlightt, int direc, double *rs,
					 double *range, double *drds);

/* mex interface -------------------------------------------------------------*/
extern void mexFunction(int nargout, mxArray *argout[], int nargin,
					    const mxArray *argin[])
{
	int corrlightt=1,direc=1;
	double *state,*rpos,*U,*dtr,*rs,*range,*drds=NULL,*pd;

	if (nargin<4||mxGetM(argin[0])<6||mxGetM(argin[1])<3||mxGetM(argin[2])<3||
		mxGetN(argin[2])<3||mxGetM(argin[3])<1)
		mexErrMsgTxt("argin error");
	if (nargout<2) mexErrMsgTxt("argout error"); 
	
	state=mxGetPr(argin[0]);
	rpos =mxGetPr(argin[1]);
	U    =mxGetPr(argin[2]);
	dtr  =mxGetPr(argin[3]);
	if (nargin>=5&&(pd=mxGetPr(argin[4]))!=NULL) corrlightt=(int)*pd;
	if (nargin>=6&&(pd=mxGetPr(argin[5]))!=NULL) direc=(int)*pd;
	
	argout[0]=mxCreateDoubleMatrix(3,1,mxREAL); rs   =mxGetPr(argout[0]);
	argout[1]=mxCreateDoubleMatrix(1,1,mxREAL); range=mxGetPr(argout[1]);
	if (nargin>2) {
		argout[2]=mxCreateDoubleMatrix(1,6,mxREAL); drds=mxGetPr(argout[2]);
	}
	SatRange(state,rpos,U,dtr,corrlightt,direc,rs,range,drds);
}
