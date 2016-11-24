/*------------------------------------------------------------------------------
% [system] : GpsTools
% [module] : satellite orbit element to position/velocity
% [func]   : transform satellite orbit elements to state(position/velocity)
%            calculate partial derivatives of state by orbit elements
% [argin]  : ele   = orbit elements(m,deg) [a,e,i,OMG,omg,M]
% [argout] : state = position/velocity(m,m/sec) [x;y;z;xdot;ydot:zdot]
%            (dsde) = partial derivatives of state by orbit elements(6x6)
% [note]   : error if a<=0,e<=0,1<=e
% [version]: $Revision: 12 $ $Date: 2008-11-25 10:02:15 +0900 (ç«, 25 11 2008) $
%            Copyright(c) 2004-2006 by T.Takasu, all rights reserved
% [history]: 04/01/05   0.1  new
*-----------------------------------------------------------------------------*/
#include "qtcmn.h"

/* matrix(3x2)xvector(2x1) ---------------------------------------------------*/
#define Mv2(X,y,z) \
{ \
	(z)[0]=(X)[0]*(y)[0]+(X)[3]*(y)[1]; \
	(z)[1]=(X)[1]*(y)[0]+(X)[4]*(y)[1]; \
	(z)[2]=(X)[2]*(y)[0]+(X)[5]*(y)[1]; \
}
/* transform satellite orbit elements to state -------------------------------*/
extern int EleToState(const double *ele, double *state, double *dsde)
{
	int k;
	double a,e,i,OMG,omg,M,E,f,rr,sf,n,r[3],v[3],PQW[9],R1[9],R2[9],R3[9];
	double ar2,ar3,drda[2],drde[2],drdM[2],dvda[2],dvde[2],dvdM[2];
	double dPQdi[6],dPQdO[6],dPQdo[6];

	a=ele[0]; e=ele[1]; i=ele[2]*M_PI/180.0; OMG=ele[3]*M_PI/180.0;
	omg=ele[4]*M_PI/180.0; M=ele[5]*M_PI/180.0;
	if (a<=0.0||e<=0.0||1.0<=e) return 0;
	
	/* solve Kepler equation */
	for (k=1,E=M;k<20;k++) {
	    f=E-e*sin(E)-M; E=E-f/(1.0-e*cos(E)); if (ABS(f)<=1E-14) break;
	}
	rr=a*(1.0-e*cos(E)); f=1.0-e*e; sf=sqrt(f); n=sqrt(GME/(a*a*a));
	r[0]=a*(cos(E)-e);
	r[1]=a*sf*sin(E);
	v[0]=-n*a*a/rr*sin(E);
	v[1]=n*a*a/rr*sf*cos(E);
	r[2]=v[2]=0.0;
	Rz(-OMG,R1); Rx(-i,R2); MM(R1,R2,R3);
	Rz(-omg,R1); MM(R3,R1,PQW);
	Mv(PQW,r,state);
	Mv(PQW,v,state+3);
	
	/* partial derivatives of state by orbit elements */
	if (dsde!=NULL) {
		ar2=a*a/(rr*rr); ar3=ar2*a/rr;
		drda[0]=r[0]/a;
		drda[1]=r[1]/a;
		drde[0]=-a-r[1]*r[1]/(rr*f);
		drde[1]=r[0]*r[1]/(rr*f);
		drdM[0]=v[0]/n;
		drdM[1]=v[1]/n;
		dvda[0]=-v[0]/(2*a);
		dvda[1]=-v[1]/(2*a);
		dvde[0]=v[0]*ar2*(2*r[0]/a+e/f*r[1]*r[1]/(a*a));
		dvde[1]=n/sf*ar2*(r[0]*r[0]/rr-r[1]*r[1]/(a*f));
		dvdM[0]=-n*ar3*r[0];
		dvdM[1]=-n*ar3*r[1];
		dPQdO[0]=-PQW[1];
		dPQdO[1]=PQW[0];
		dPQdO[2]=dPQdO[5]=0.0;
		dPQdO[3]=-PQW[4];
		dPQdO[4]=PQW[3];
		for (k=0;k<3;k++) {
			dPQdi[k]=sin(omg)*PQW[6+k];
			dPQdi[3+k]=cos(omg)*PQW[6+k];
			dPQdo[k]=PQW[3+k];
			dPQdo[3+k]=-PQW[k];
		}
		Mv2(PQW,drda,dsde);		/* ds/da */
		Mv2(PQW,dvda,dsde+3);
		Mv2(PQW,drde,dsde+6);	/* ds/de */
		Mv2(PQW,dvde,dsde+9);
		Mv2(dPQdi,r,dsde+12);	/* ds/di */
		Mv2(dPQdi,v,dsde+15);
		Mv2(dPQdO,r,dsde+18);	/* ds/dƒ¶ */
		Mv2(dPQdO,v,dsde+21);
		Mv2(dPQdo,r,dsde+24);	/* ds/dƒÖ */
		Mv2(dPQdo,v,dsde+27);
		Mv2(PQW,drdM,dsde+30);	/* ds/dM */
		Mv2(PQW,dvdM,dsde+33);
	}
	return 1;
}
