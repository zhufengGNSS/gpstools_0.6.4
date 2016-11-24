/*------------------------------------------------------------------------------
% [system] : GpsTools
% [module] : navigation messages to satellite position/velocity
% [func]   : calculate satllite position/velocity from navigation messages
% [argin]  : td    = date (mjd-gpst)
%            ts    = time (sec)
%            nav   = navigation messages
%           (type) = GPS/GLONASS('N'=GPS,'G'=GLONASS)
%           (opt)  = option (1:relativity correction off)
% [argout] : pos = satellite position [x;y;z] (m)(ecef)
%            dts = satellite clock error [bias;drift;drift-rate]
%                  bias/drift/drift-rate(sec,sec/sec,sec/sec^2)
%           (vel)= satellite velocity [vx;vy;vz] (m/sec)(ecef)
%           (svh)= satellite health
% [note]   :
% [version]: $Revision: 16 $ $Date: 2008-12-12 15:49:30 +0900 (金, 12 12 2008) $
%            Copyright(c) 2004-2008 by T.Takasu, all rights reserved
% [history]: 04/11/16  0.1  new
%            06/02/28  0.2  add argin opt
%            08/12/04  0.3  add argout svh (gt_0.6.4)
%                           suppress warning
%-----------------------------------------------------------------------------*/
#include "mex.h"

extern int NavToState(double td, double ts, const double *nav, int nnav,
				      char type, int opt, double *pos, double *dts, double *vel);

/* mex interface -------------------------------------------------------------*/
extern void mexFunction(int nargout, mxArray *argout[], int nargin,
					    const mxArray *argin[])
{
	double *td,*ts,*nav,*pos,*dts,*vel=NULL;
	char type[2]="N";
	int i,nnav,opt=0,svh=0;
	
	if (nargin<3||mxGetM(argin[0])<1||mxGetM(argin[1])<1)
		mexErrMsgTxt("argin error");
	
	if (nargout>4) mexErrMsgTxt("argout error"); 
	
	td =mxGetPr(argin[0]);
	ts =mxGetPr(argin[1]);
	nav=mxGetPr(argin[2]); nnav=(int)mxGetM(argin[2]);
	if (nnav>0&&mxGetN(argin[2])!=37) mexErrMsgTxt("argin error");
	if (nargin>=4&&mxIsChar(argin[3]))
		mxGetString(argin[3],type,sizeof(type));
	if (nargin>=5&&mxGetM(argin[4])>=1) opt=(int)*mxGetPr(argin[4]);
	
	argout[0]=mxCreateDoubleMatrix(3,1,mxREAL); pos=mxGetPr(argout[0]);
	argout[1]=mxCreateDoubleMatrix(3,1,mxREAL); dts=mxGetPr(argout[1]);
	if (nargout>=3) {
		argout[2]=mxCreateDoubleMatrix(3,1,mxREAL);
		vel=mxGetPr(argout[2]);
	}
	if (nnav>0) {
		svh=NavToState(*td,*ts,nav,nnav,type[0],opt,pos,dts,vel);
	}
	else {
		for (i=0;i<3;i++) {
			pos[i]=mxGetNaN();
			dts[i]=mxGetNaN();
			if (vel) vel[i]=mxGetNaN();
		}
	}
	if (nargout>=4) {
		argout[3]=mxCreateDoubleMatrix(1,1,mxREAL);
		*mxGetPr(argout[3])=(double)svh;
	}
}
