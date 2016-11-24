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
%           (dts)= satellite clock error [bias;drift;drift-rate]
%                  bias/drift/drift-rate(sec,sec/sec,sec/sec^2)
%           (vel)= satellite velocity [vx;vy;vz] (m/sec)(ecef)
% [note]   :
% [version]: $Revision: 5 $ $Date: 06/07/08 1:01 $
%            Copyright(c) 2004-2006 by T.Takasu, all rights reserved
% [history]: 04/11/16  0.1  new
%            06/02/28  0.2  add argin opt
%-----------------------------------------------------------------------------*/
#include "mex.h"
#include "qtcmn.h"

extern void CalToMjd(const double *dt, double *mjd, double *sec, int len_dt);

/* gps navigation message to satellite postion -------------------------------*/
static void GpsNavToPos(const double *nav, double tk, double *pos, double *ek)
{
	static const double mu=3.986005E14;
	static const double omge=7.2921151467E-5;
	double a,e,n,omg,mk,ekk,phik,dk[3],uk,rk,ik,omgk;
	double cos2p,sin2p,cosu,sinu,coso,sino;
	int i;
	if (tk<-302400.0) tk=tk+604800.0; else if (tk>302400.0) tk=tk-604800.0;
	a=nav[16]*nav[16]; e=nav[14]; n=sqrt(mu/(a*a*a))+nav[11]; omg=nav[23];
	mk=nav[12]+n*tk;
	for (i=0,*ek=mk;i<30;i++) {
		ekk=mk+e*sin(*ek); if (ABS(ekk-*ek)<1E-12) break; else *ek=ekk;
	}
	phik=atan2(sqrt(1-e*e)*sin(*ek),cos(*ek)-e)+omg;
	cos2p=cos(2*phik); sin2p=sin(2*phik);
	dk[0]=nav[13]*cos2p+nav[15]*sin2p;
	dk[1]=nav[22]*cos2p+nav[10]*sin2p;
	dk[2]=nav[18]*cos2p+nav[20]*sin2p;

	uk=phik+dk[0];
	rk=a*(1-e*cos(*ek))+dk[1];
	ik=nav[21]+nav[25]*tk+dk[2];
	omgk=nav[19]+(nav[24]-omge)*tk-omge*nav[17];
    cosu=cos(uk); sinu=sin(uk); coso=cos(omgk); sino=sin(omgk);
	pos[0]=rk*(coso*cosu-cos(ik)*sino*sinu);
	pos[1]=rk*(sino*cosu+cos(ik)*coso*sinu);
	pos[2]=rk*sin(ik)*sinu;
}
/* gps navigation message to satellite position/clock/velocity ---------------*/
static void GpsNavToState(const double *nav, double t_toe, double t_toc,
						  int opt, double *pos, double *dts, double *vel)
{
	double posm[3],posp[3],dt=1E-3,ek;
	
	/* satellite position */
	GpsNavToPos(nav,t_toe,pos,&ek);

	/* satellite clock */
	if (dts!=NULL) {
		dts[0]=nav[6]+nav[7]*t_toc+nav[8]*t_toc*t_toc/2.0;
		dts[1]=nav[7]+nav[8]*t_toc;
		dts[2]=nav[8];
		
		/* relativistic correction */
		if (!opt) dts[0]-=4.443E-10*nav[14]*nav[16]*sin(ek);
	}
	/* satellite velocity by differential approximation */
	if (vel!=NULL) {
		GpsNavToPos(nav,t_toe-dt,posm,&ek);
		GpsNavToPos(nav,t_toe+dt,posp,&ek);
		vel[0]=(posp[0]-posm[0])/(2.0*dt);
		vel[1]=(posp[1]-posm[1])/(2.0*dt);
		vel[2]=(posp[2]-posm[2])/(2.0*dt);
	}
}
/* gps navigation messages to satellite position/clock/velocity --------------*/
extern void NavToState(double td, double ts, const double *nav, int nnav,
					   char type, int opt, double *pos, double *dts, double *vel)
{
	double t_toe,t_toc,t0,dt,navi[37],dtmin=99999999.0;
	int i,imin=0;
	if (type=='N') {
		/* search navigation message having closest toe to the time */
		for (i=imin=0;i<nnav;i++) {
			dt=(td-44244.0-nav[i+27*nnav]*7)*86400.0+ts-nav[i+17*nnav];
			if (ABS(dt)<ABS(dtmin)) {dtmin=dt; imin=i;}
		}
		t_toe=dtmin;
		for (i=0;i<37;i++) navi[i]=nav[imin+i*nnav];
		navi[0]+=nav[0]<70.0?2000.0:1900.0;
		CalToMjd(navi,&t0,NULL,3);
		t_toc=(td-t0)*86400.0+ts-navi[3]*3600.0-navi[4]*60.0+navi[5];
		
		/* navigation message to satellite position/clock/velocity */
		GpsNavToState(navi,t_toe,t_toc,opt,pos,dts,vel);
	}
	else if (type=='G') {
		mexErrMsgTxt("GLONASS not supported yet\n");
	}
}
