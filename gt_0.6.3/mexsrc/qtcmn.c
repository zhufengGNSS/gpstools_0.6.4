/*------------------------------------------------------------------------------
% [system] : GpsTools
% [module] : common routines
% [func]   : common routines
% [argin]  : none
% [argout] : none
% [note]   : without external library
% [version]: $Revision: 5 $ $Date: 06/07/08 1:01 $
%            Copyright(c) 2004-2006 by T.Takasu, all rights reserved
% [history]: 04/01/05   0.1  new
%-----------------------------------------------------------------------------*/
#include "qtcmn.h"

/* product of matrixes -------------------------------------------------------*/
extern void MatMul(const double *X, const double *Y, int m, int n, int k,
				   double *Z)
{
	int i,j,p;
	const double *x,*y;
	double *z,w;

	memset(Z,0,sizeof(double)*m*k);

	for (i=0,y=Y,z=Z;i<k;i++,y+=n,z+=m) {
		for (j=0,x=X;j<n;j++,x+=m) {
			for (p=0,w=y[j];p<m;p++) {
				z[p]+=x[p]*w;
			}
		}
	}
}
/* LU decomposition ----------------------------------------------------------*/
static int LU(const double *X, int n, double *Y, int *index)
{
	int i,j,k,imax=0;
	const double *x;
	double *v,w,max,*y,*yy;

	if ((v=(double *)malloc(sizeof(double)*n))==NULL) return 0;
	for (i=0,y=Y;i<n;i++,y+=n) for (j=0,x=X+i;j<n;j++,x+=n) y[j]=*x;
	
	for (i=0,y=Y;i<n;i++,y+=n) {
		for (j=0,max=0.0;j<n;j++) max=MAX(ABS(y[j]),max);
		if (max>0.0) v[i]=1.0/max; else {free(v); return 0;}
	}
	for (j=0;j<n;j++) {
		for (i=0,y=Y;i<j;i++,y+=n)
			for (k=0,yy=Y+j;k<i;k++,yy+=n) y[j]-=y[k]*(*yy);
		
		for (i=j,y=Y+i*n,max=0.0;i<n;i++,y+=n) {
			w=y[j];
			for (k=0,yy=Y+j;k<j;k++,yy+=n) w-=y[k]*(*yy);
			y[j]=w;
			if ((w=v[i]*ABS(y[j]))>=max) {max=w; imax=i;}
		}
		y=Y+j*n;
		if (j!=imax) {
			for (k=0,yy=Y+imax*n;k<n;k++) SWAP(y[k],yy[k]);
			v[imax]=v[j];
		}
		index[j]=imax;
		if (ABS(y[j])<=0.0) {free(v); return 0;}
		if (j!=n-1) for (i=j+1,yy=Y+i*n+j;i<n;i++,yy+=n) *yy/=y[j];
	}
	free(v);
	return 1;
}
/* inverse matrix ------------------------------------------------------------*/
extern int MatInv(const double *X, int n, double *Y)
{
	int i,j,k,jj,*index;
	double *D,*x,*y,w;

	D=(double *)malloc(sizeof(double)*n*n);
	index=(int *)malloc(sizeof(int)*n);
	if (!D||!index) {free(D); free(index); return 0;}

	if (!LU(X,n,D,index)) {free(D); free(index); return 0;}
	
	for (i=0,y=Y;i<n;i++,y+=n) {
		memset(y,0,sizeof(double)*n); y[i]=1.0;
		for (j=0,x=D,jj=-1;j<n;j++,x+=n) {
			SWAP(y[index[j]],y[j]);
			if (jj>=0) {
				w=y[j];
				for (k=0;k<j;k++) w-=x[k]*y[k];
				y[j]=w;
			}
			else if (y[j]!=0.0) jj=j;
		}
		for (j--,x-=n;j>=0;j--,x-=n) {
			w=y[j];
			for (k=n-1;k>j;k--) w-=x[k]*y[k];
			y[j]=w/x[j];
		}
	}
	free(D); free(index);
	return 1;
}
/* transpose matrix ----------------------------------------------------------*/
extern void MatTr(const double *X, int m, int n, double *Y)
{
	int i,j;
	const double *x;
	double *y;
	for (i=0,x=X;i<n;i++) {
		for (j=0,y=Y+i;j<m;j++,x++,y+=n) {
			*y=*x;
		}
	}
}
