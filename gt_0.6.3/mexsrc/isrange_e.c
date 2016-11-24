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
#include "qtcmn.h"

/* satellite postion approximation(t+dt) -------------------------------------*/
static void SatPos(const double *state, double dt, double *rs)
{
	double r=NORM(state),k=1.0-GME/(2.0*r*r*r)*dt*dt;
	rs[0]=state[0]*k+state[3]*dt;
	rs[1]=state[1]*k+state[4]*dt;
	rs[2]=state[2]*k+state[5]*dt;
}
/* transition matrix of position(t to t+dt) ----------------------------------*/
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
/* range between satellites --------------------------------------------------*/
extern void IsRange(const double *state1, const double *state2,
				    const double *dts1, const double *dts2, int corrlightt,
				    double *rs1, double *rs2, double *range, double *drds1,
				    double *drds2)
{
	int i,j;
	double dt,rng,rss[3],phirr[9],R1[9];
	
	if (!corrlightt) {
		for (i=0;i<3;i++) {
			rs1[i]=state1[i];
			rs2[i]=state2[i];
			rss[i]=rs1[i]-rs2[i];
		}
		range[0]=NORM(rss);
		if (drds1!=NULL) {
			for (i=0;i<3;i++) drds1[i]=rss[i]/range[0];
			for (i=3;i<6;i++) drds1[i]=0.0;
		}
		if (drds2!=NULL) {
			for (i=0;i<3;i++) drds2[i]=-rss[i]/range[0];
			for (i=3;i<6;i++) drds2[i]=0.0;
		}
	}
	else { /* use light time/receiver clock bias corrections */
		SatPos(state2,-dts2[0],rs2);
		range[0]=0.0;
		for (i=0;i<10;i++) {
			dt=dts2[0]+range[0]/M_C;
			SatPos(state1,-dt,rs1);
			for (j=0;j<3;j++) rss[j]=rs1[j]-rs2[j];
			rng=range[0]; range[0]=NORM(rss);
			if (ABS(range[0]-rng)<1E-4) break;
		}
		for (i=0;i<3;i++) rss[i]/=range[0];
		if (drds1!=NULL) {
			Phirr(state1,-dt,phirr);
			Tr(phirr,R1);
			Mv(R1,rss,drds1);
			for (i=0;i<3;i++) drds1[3+i]=-rss[i]*dt;
		}
		if (drds2!=NULL) {
			Phirr(state2,-dts2[0],phirr);
			Tr(phirr,R1);
			for (i=0;i<3;i++) rss[i]=-rss[i];
			Mv(R1,rss,drds2);
			for (i=0;i<3;i++) drds2[3+i]=rss[i]*dts2[0];
		}
	}
}
