/*------------------------------------------------------------------------------
% [system] : GpsTools
% [module] : satellite state to orbit element
% [func]   : convert satellite position/velocity to orbit element
% [argin]  : state = position/velocity(m,m/sec) [x;y;z;xdot;ydot:zdot]
% [argout] : ele   = orbit element(m,deg) [a,e,i,OMG,omg,M]
% [note]   : tangental orbit element
% [version]: $Revision: 3 $ $Date: 06/07/08 1:01 $
%            Copyright(c) 2004-2006 by T.Takasu, all rights reserved
% [history]: 04/01/05   0.1  new
%----------------------------------------------------------------------------*/
#include "mex.h"
#include "qtcmn.h"

/* satellite state to orbit element -----------------------------------------*/
extern void StateToEle(const double *state, double *ele)
{
	int j;
	double r[3],v[3],h[3],w[3],hh,rr,aa,a,e,i,OMG,omg,M,E;
	
	r[0]=state[0]; r[1]=state[1]; r[2]=state[2];
	v[0]=state[3]; v[1]=state[4]; v[2]=state[5];
	rr=NORM(r);
	CROSS(r,v,h);
	aa=2.0/rr-DOT(v,v)/GME;
	if (NORM(h)<1E-12||aa<=0) {
		for (j=0;j<6;j++) ele[j]=mxGetNaN(); /* no rotating orbit */
		return;
	}
	a=1.0/aa;
	hh=NORM(h);
	w[0]=h[0]/hh; w[1]=h[1]/hh; w[2]=h[2]/hh;
	i=atan(NORM2(w)/w[2]);
	OMG=atan2(h[0],-h[1]);
	e=sqrt(1.0-DOT(h,h)/(GME*a));
	E=atan2(DOT(r,v)/sqrt(GME*a),1.0-rr/a);
	M=E-e*sin(E);
	omg=atan2(r[2],-r[0]*w[1]+r[1]*w[0])-atan2(sqrt(1-e*e)*sin(E),cos(E)-e);
	if (i<0.0) i=i+M_PI;
	if (OMG<0.0) OMG=OMG+2.0*M_PI;
	if (omg<0.0) omg=omg+2.0*M_PI;
	if (M<0.0) M=M+2.0*M_PI;
	ele[0]=a;
	ele[1]=e;
	ele[2]=i*180.0/M_PI;
	ele[3]=OMG*180.0/M_PI;
	ele[4]=omg*180.0/M_PI;
	ele[5]=M*180.0/M_PI;
}
