/*-------------------------------------------------------------------------------
% [system] : GpsTools
% [module] : range between satellites
% [func]   : calculate range between two satellites
% [argin]  : state1 = satellite no.1 position/velocity(m) (eci)
%            state2 = satellite no.2 position/velocity(m) (eci)
%            dts1   = satellite no.1 clock bias(sec)
%            dts2   = satellite no.2 clock bias(sec)
%            corrlightt = light-time/receiver clock bias correction flag (1:on)
% [argout] : rs1    = satellite no.1 tx position (m) (eci)
%            rs2    = satellite no.2 tx position (m) (eci)
%            range  = range between satellite no.1 and no.2
%            (drds1)= partial derivatives of range by satellite no.1 state
%            (drds2)= partial derivatives of range by satellite no.2 state
% [note]   :
% [version]: $Revision: 4 $ $Date: 06/07/08 1:01 $
%            Copyright(c) 2004-2006 by T.Takasu, all rights reserved
% [history]: 04/02/11   0.1  new
%-----------------------------------------------------------------------------*/
#include "mex.h"

extern void IsRange(const double *state1, const double *state2,
				    const double *dts1, const double *dts2, int corrlightt,
				    double *rs1, double *rs2, double *range, double *drds1,
				    double *drds2);

/* mex interface -------------------------------------------------------------*/
extern void mexFunction(int nargout, mxArray *argout[], int nargin,
					    const mxArray *argin[])
{
	int corrlightt;
	double *state1,*state2,*dts1,*dts2,*rs1,*rs2,*range,*drds1=NULL,*drds2=NULL;

	if (nargin<5||mxGetM(argin[0])<6||mxGetM(argin[1])<6||mxGetM(argin[2])<1||
		mxGetM(argin[3])<1||mxGetM(argin[4])<1)
		mexErrMsgTxt("argin error");
	if (nargout<3) mexErrMsgTxt("argout error"); 
	
	state1=mxGetPr(argin[0]);
	state2=mxGetPr(argin[1]);
	dts1  =mxGetPr(argin[2]);
	dts2  =mxGetPr(argin[3]);
	corrlightt=(int)*mxGetPr(argin[4]);
	
	argout[0]=mxCreateDoubleMatrix(3,1,mxREAL); rs1=mxGetPr(argout[0]);
	argout[1]=mxCreateDoubleMatrix(3,1,mxREAL); rs2=mxGetPr(argout[1]);
	argout[2]=mxCreateDoubleMatrix(1,1,mxREAL); range=mxGetPr(argout[2]);
	if (nargin>3) {
		argout[3]=mxCreateDoubleMatrix(1,6,mxREAL); drds1=mxGetPr(argout[3]);
	}
	if (nargin>4) {
		argout[4]=mxCreateDoubleMatrix(1,6,mxREAL); drds2=mxGetPr(argout[4]);
	}
	IsRange(state1,state2,dts1,dts2,corrlightt,rs1,rs2,range,drds1,drds2);
}
