/*------------------------------------------------------------------------------
% [system] : GpsTools
% [module] : double difference
% [func]   : generate double difference between satellites/stations
% [argin]  : dz  = residuals (O-C) (m)
%            dgds= partial derivatives
%            sig = measurement noise std. deviation (m)
%            ig  = observation index [t,isat,ircv,arcf;...]
%            ircv= station numbers of baselines [ircv1,ircv2;...]
% [argout] : dzs = double diff. of resuduals
%            G   = double diff. of partial derivatives
%            R   = covariance matrix of measurement noise
%            iz  = satellite/station index[isat1,isat2,ircv1,ircv2;...]
% [note]   :
% [version]: $Revision: 3 $ $Date: 06/07/08 1:01 $
%            Copyright(c) 2004-2006 by T.Takasu, all rights reserved
% [history]: 04/06/19  0.1  new
%            05/12/02  0.2  select highest elevation satellite as reference
%-----------------------------------------------------------------------------*/
#include "mex.h"
#include "qtcmn.h"

#define		NSATMAX		32		/* max satellite number */
#define		NZMAX		1024	/* max observation data */
#define		NZSMAX		512		/* max double difference */

/* double difference ---------------------------------------------------------*/
extern int DblDiff(const double *dz, const double *dgds, const double *sig,
          	       const double *ig, const double *ircv, int nz, int nx, int nr,
          	       double *dzs, double *G, double *R, double *iz)
{
	const double *is=ig+nz,*ir=ig+2*nz;
	int i,j,k,j1,j2,nzs,iz1[NSATMAX],iz2[NSATMAX],id[4*NZSMAX],*p=id,w[NZMAX];
	double var;
	
	for (i=nzs=0;i<nr;i++) { /* each station pairs */
		for (j=0;j<NSATMAX;j++) iz1[j]=iz2[j]=-1;
		for (j=0;j<nz;j++) {
			if ((int)ir[j]==(int)ircv[i]) iz1[(int)is[j]]=j;
			else if ((int)ir[j]==(int)ircv[i+nr]) iz2[(int)is[j]]=j;
		}
		for (j1=-1,j2=0;j2<NSATMAX;j2++) {
			
			if (iz1[j2]<0||iz2[j2]<0) continue;
			
			if (j1>=0) { /* each satellite pairs */
				p[0]=iz1[j1];
				p[1]=iz1[j2];
				p[2]=iz2[j1];
				p[3]=iz2[j2];
				dzs[nzs]=(dz[p[0]]-dz[p[1]])-(dz[p[2]]-dz[p[3]]); /* double diff. */
				for (j=0;j<nx;j++)
					G[nzs+j*NZSMAX]=(dgds[p[0]+j*nz]-dgds[p[1]+j*nz])-
									(dgds[p[2]+j*nz]-dgds[p[3]+j*nz]);
				iz[nzs]=(double)j1;
				iz[nzs+NZSMAX]=(double)j2;
				iz[nzs+NZSMAX*2]=ircv[i];
				iz[nzs+NZSMAX*3]=ircv[i+nr];
				nzs++; p+=4;
				if (nzs>=NZSMAX) mexErrMsgTxt("dbldiff size overflow\n");
			}
			j1=j2;
		}
	}
	for (i=0;i<nzs;i++) /* cov. matrix of measurement noise */
	for (j=0;j<=i;j++) {
		for (k=0;k<nz;k++) w[k]=0;
		for (k=0;k<4;k++) {w[id[i*4+k]]+=1; w[id[j*4+k]]+=1;}
		for (k=0,var=0.0;k<nz;k++) if (w[k]==2) var+=sig[k]*sig[k];
		R[i*NZSMAX+j]=R[j*NZSMAX+i]=var;
	}
	return nzs;
}
