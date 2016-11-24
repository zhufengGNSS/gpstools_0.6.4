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
%            (drds) = partial derivatives of range by satellite state
% [note]   :
% [version]: $Revision: 12 $ $Date: 2008-11-25 10:02:15 +0900 (火, 25 11 2008) $
%            Copyright(c) 2004-2006 by T.Takasu, all rights reserved
% [history]: 04/04/12   0.1  new
%-----------------------------------------------------------------------------*/
#include "qtcmn.h"

/* satellite position(t+dt) --------------------------------------------------*/
static void SatPos(const double *state, double dt, double *rs)
{
	double r=NORM(state),k=1.0-GME/(2.0*r*r*r)*dt*dt;
	rs[0]=state[0]*k+state[3]*dt;
	rs[1]=state[1]*k+state[4]*dt;
	rs[2]=state[2]*k+state[5]*dt;
}
/* station postion(t+dt) -----------------------------------------------------*/
#define RcvPos(rpos,U,dt,rr) \
{ \
	double _R1[9],_R2[9],_R3[9]; \
	Tr(U,_R1); Rz(-(dt)*0.7292115E-4,_R2); MM(_R1,_R2,_R3); Mv(_R3,rpos,rr); \
}
/* state transition matrix(t to t+dt) ----------------------------------------*/
static void Phirr(const double *state, double dt, double *phirr)
{
	double r=NORM(state),r2=r*r,k=GME/(2.0*r2*r2*r)*dt*dt;
	phirr[0]=1.0+k*(3.0*state[0]*state[0]-r2);
	phirr[4]=1.0+k*(3.0*state[1]*state[1]-r2);
	phirr[8]=1.0+k*(3.0*state[2]*state[2]-r2);
	phirr[1]=phirr[3]=k*3.0*state[0]*state[1];
	phirr[2]=phirr[6]=k*3.0*state[0]*state[2];
	phirr[5]=phirr[7]=k*3.0*state[1]*state[2];
}
/* satellite-station range ---------------------------------------------------*/
extern void SatRange(const double *state, const double *rpos, const double *U,
					 const double *dtr, int corrlightt, int direc, double *rs,
					 double *range, double *drds)
{
	int i,j;
	double rng,dt,urpos[3],rsr[3],rr[3],phirr[9],R1[9];
	
	/* no light-time correction */
	if (!corrlightt) {
		for (i=0;i<3;i++) rs[i]=state[i];
		Tr(U,R1); Mv(R1,rpos,urpos);
		for (i=0;i<3;i++) rsr[i]=rs[i]-urpos[i];
		range[0]=NORM(rsr);
		if (drds!=NULL) {
			for (i=0;i<3;i++) drds[i]=rsr[i]/range[0];
			for (i=3;i<6;i++) drds[i]=0.0;
		}
	}
	/* station to satellite */
	else if (direc==2) {
		SatPos(state,-dtr[0],rs);
		range[0]=0.0;
	    for (i=0;i<10;i++) {
	        dt=dtr[0]+range[0]/M_C;
	        RcvPos(rpos,U,-dt,rr);
	        for (j=0;j<3;j++) rsr[j]=rr[j]-rs[j];
	        rng=range[0]; range[0]=NORM(rsr);
	        if (ABS(range[0]-rng)<1E-4) break;
	    }
		if (drds!=NULL) {
			for (i=0;i<3;i++) rsr[i]/=range[0];
			Phirr(state,-dtr[0],phirr);
			Tr(phirr,R1);
			Mv(R1,rsr,drds);
			for (i=0;i<3;i++) drds[3+i]=-rsr[i]*dtr[0];
		}
	}
	/* satellite to station */
	else {
		RcvPos(rpos,U,-dtr[0],rr);
		range[0]=0.0;
	    for (i=0;i<10;i++) {
	        dt=dtr[0]+range[0]/M_C;
	        SatPos(state,-dt,rs);
	        for (j=0;j<3;j++) rsr[j]=rs[j]-rr[j];
	        rng=range[0]; range[0]=NORM(rsr);
	        if (ABS(range[0]-rng)<1E-4) break;
	    }
		if (drds!=NULL) {
			for (i=0;i<3;i++) rsr[i]/=range[0];
			Phirr(state,-dt,phirr);
			Tr(phirr,R1);
			Mv(R1,rsr,drds);
			for (i=0;i<3;i++) drds[3+i]=-rsr[i]*dt;
		}
	}
}
