/*-------------------------------------------------------------------------------
% [system] : GpsTools
% [module] : point positioning
% [argin]  : td   = date(mjd-gpst)
%            ts   = time(sec)
%            z    = observation data(pseudo-range)
%            iz   = observation data satellite/station index
%            nav  = navigation messages
%            inav = navigation messages satellite index
%            posr = approx station position (m) (ecef)
%            cdtr = approx receiver clock bias (m)
%           (dflg)= debug flag (1:on,0:off)
% [argout] : posr = station postion (m) (ecef)
%            cdtr = receiver clock bias (m)
% [note]   :
% [version]: $Revision: 13 $ $Date: 06/07/08 14:02 $
%            Copyright(c) 2004-2006 by T.Takasu, all rights reserved
% [history]: 04/10/19  0.1  new
%            05/06/07  0.2  prevent to switch ephemeris msg. in iteration
%            05/07/13  0.3  change outlier exclusion algorithm
%            05/08/08  0.4  add least square adjustment iteration
%            05/08/14  0.5  add argin dflg, add earth rotation effect
%            05/09/16  0.6  return nan if obs data count < 4
%-----------------------------------------------------------------------------*/
#include "qtcmn.h"
#include "mex.h"

#define	MINOBS		4				/* min observation data count */
#define	MAXRMS		10.0			/* max residuals rms (m) */
#define	OMGE		0.7292115E-4	/* earth rotation rate (rad/sec) */

extern void EcefToGeod(const double *epos, double *gpos, double *E);
extern void NavToState(double td, double ts, const double *nav, int nnav,
					   char type, int opt, double *pos, double *dts, double *vel);
extern void SatAzEl(const double *spos, const double *rpos, double *azel);
extern void trop_saast(const double *t, const double *azel, const double *gpos,
					   const double *met_prm, double *trop);

/* least square estimation ---------------------------------------------------*/
static void Lsq(const double *dz, const double *H, int nx, int nz, double *dx)
{
	double *Ht,*HH,*Hdz,*invHH;
	Ht=MAT(nz,nx); HH=MAT(nx,nx); Hdz=VEC(nx); invHH=MAT(nx,nx);
	if (!Ht||!HH||!Hdz||!invHH) mexErrMsgTxt("malloc error\n");
	MatTr(H,nx,nz,Ht);
	MatMul(H,Ht,nx,nz,nx,HH);
	MatMul(H,dz,nx,nz,1,Hdz);
	MatInv(HH,nx,invHH);
	MatMul(invHH,Hdz,nx,nx,1,dx);
	FreeMat(Ht); FreeMat(HH); FreeMat(Hdz); FreeMat(invHH);
}
/* calculate residuals -------------------------------------------------------*/
static int residual(const double *td, const double *ts, const double *z,
					const double *iz, int *vz, int nz, const double *nav,
					const double *inav, int nnav, const double *posr,
					const double *cdtr, const double *dflg, double *dz, double *H)
{
	double gpos[3],poss[3],posz[3],dts[3],r[3],azel[2],navi[37],R[9];
	double range,rng,trop,dt,dtmin;
	int i,j,imin,n=0;
	
	EcefToGeod(posr,gpos,NULL);
	
	for (i=0;i<nz;i++) {
		if (!vz[i]) continue;
		
		/* search navigation message having closest toe to the time */
		for (j=0,imin=-1,dtmin=99999999.0;j<nnav;j++) {
			if ((int)inav[j]!=(int)iz[i]) continue;
			dt=(td[0]-44244.0-nav[j+27*nnav]*7)*86400.0+ts[0]-nav[j+17*nnav];
			if (ABS(dt)<ABS(dtmin)) {dtmin=dt; imin=j;}
		}
		if (imin==-1) {
		    /*mexPrintf("warning : no navigation msg : ts=%.0f,sat=%2.0f\n",*ts,iz[i]);*/
			vz[i]=0;
			continue;
		}
		for (j=0;j<37;j++) navi[j]=nav[imin+j*nnav];

		/* light-time iteration */
		for (j=0,rng=0.0;j<10;j++) {
			/* satellite position */
			NavToState(td[0],ts[0]-(cdtr[0]+rng)/M_C,navi,1,'N',0,poss,dts,NULL);
			
			/* earth rotation effect */
			Rz(OMGE*(cdtr[0]+rng)/M_C,R);
			Mv(R,poss,posz);
			
			/* satellite range */
			r[0]=posr[0]-posz[0];
			r[1]=posr[1]-posz[1];
			r[2]=posr[2]-posz[2];
			range=NORM(r);
			if (ABS(range-rng)<1E-4) break;
			rng=range;
		}
		/* satellite azimath/elevation angle */
		SatAzEl(poss,posr,azel);
		
		/* tropospheric delay */
		trop_saast(td,azel,gpos,NULL,&trop);
		
		/* residual */
		dz[n]=z[i]-(range+cdtr[0]-M_C*dts[0]+trop);
		
		/* partial derivatives */
		H[n*4  ]=r[0]/range;
		H[n*4+1]=r[1]/range;
		H[n*4+2]=r[2]/range;
		H[n*4+3]=1.0;
		if ((int)*dflg>=2)
			mexPrintf("t=%6.0f,sat=%2.0f,dz=%.3f,z=%.3f,h=%.3f,%.3f,%.3f,%.3f\n",
					  *ts,iz[i],dz[n],z[i],range,cdtr[0],M_C*dts[0],trop);
		n++;
	}
	return n;
}
/* least square adjustment and iteration -------------------------------------*/
static double AdjustPos(const double *td, const double *ts, const double *z,
						const double *iz, int *vz, int nz, const double *nav,
						const double *inav, int nnav, const double *dflg,
						double *x)
{
	double dx[4],*H,*dz,rms;
	int i,j,n;

	H=MAT(4,nz); dz=VEC(nz); if (!H||!dz) mexErrMsgTxt("malloc error\n");
	
	for (i=0;i<10;i++) {
		/* residulas */
		if ((n=residual(td,ts,z,iz,vz,nz,nav,inav,nnav,x,x+3,dflg,dz,H))<MINOBS)
			break;
		
		/* least square adjustment */
		Lsq(dz,H,4,n,dx); for (j=0;j<4;j++) x[j]+=dx[j];
		
		if (ABS(dx[0])<0.1&&ABS(dx[1])<0.1&&ABS(dx[2])<0.1&&ABS(dx[3])<0.1)
			break;
	}
	if (i>=10) rms=9999999.9; /* divergent */
	else if (n<MINOBS) rms=-1.0; /* lack of obs data count */
	else {
		/* residual rms */
		for (i=0,rms=0.0;i<n;i++) rms+=dz[i]*dz[i];
		rms=sqrt(rms/(double)n);
	}
	FreeMat(H); FreeMat(dz);
	if ((int)*dflg>=2) mexPrintf("t=%6.0f,rms=%.3f\n",*ts,rms);
	
	return rms;
}
/* point positioning ---------------------------------------------------------*/
extern void PointPos(const double *td, const double *ts, const double *z,
					 const double *iz, int nz, const double *nav,
					 const double *inav, int nnav, const double *posr0,
					 const double *cdtr0, const double *dflg, double *posr,
					 double *cdtr)
{
	int i,j,n,vz[32];
	double x[4],rms,r,minr;
	
	for (i=0;i<nz;i++) vz[i]=1;
	for (i=0;i<3;i++) posr[i]=mxGetNaN(); *cdtr=mxGetNaN();
	if (nz<=0||nnav<=0) return;
	
	while (1) {
		/* least square adjustment */
		for (i=0;i<3;i++) x[i]=posr0[i]; x[3]=*cdtr0;
		if ((rms=AdjustPos(td,ts,z,iz,vz,nz,nav,inav,nnav,dflg,x))<0.0) return;
		if (rms<MAXRMS) break;
		
		/* exclude outlier */
		for (i=n=0,minr=999999999.9;i<nz;i++) {
			if (!vz[i]) continue;
			vz[i]=0;
			for (j=0;j<3;j++) x[j]=posr0[j]; x[3]=*cdtr0;
			if ((r=AdjustPos(td,ts,z,iz,vz,nz,nav,inav,nnav,dflg,x))<minr) {
				minr=r; n=i;
			}
			vz[i]=1;
		}
		vz[n]=0;
		if ((int)*dflg>=1)
			mexPrintf("t=%6.0f,rms=%.3f excluded sat=%2.0f\n",*ts,rms,iz[n]);
	}
	if ((int)*dflg>=1) mexPrintf("t=%6.0f,rms=%.3f\n",*ts,rms);
	
	for (i=0;i<3;i++) posr[i]=x[i]; *cdtr=x[3];
}
