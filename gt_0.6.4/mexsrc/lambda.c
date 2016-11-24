/*------------------------------------------------------------------------------
% [system] : GpsTools
% [module] : LAMBDA/MLAMBDA integer least-square estimation
% [func]   : LAMBDA/MLAMBDA integer least-square estimation
% [argin]  : a = float (real) solution
%            Q = covariance of float solution
%            n = number of fixed solutions
%            (flg) = reduction flag (0:off,1:on) (default:on)
% [argout] : a = fixed (integer) solutions [a1,a2,..,an]
%            f = distances of fixed solutions [f1,f2,...,fn] (f1<f2<...<fn)
%           (Z)= Z-transformation matrix
% [note]   : reference
%            X.W.Chang et al., MLAMBDA: A modified LAMBDA method for integer
%            least-squares estimation, ION AM, 2005
% [version]: $Revision: 12 $ $Date: 2008-11-25 10:02:15 +0900 (火, 25 11 2008) $
%            Copyright(c) 2004-2006 by T.Takasu, all rights reserved
% [history]: 05/12/13  0.1  new
%            06/01/27  0.2  add argin flg
%-----------------------------------------------------------------------------*/
#include "mkl_types.h"
#include "mkl_blas.h"
#include "mkl_lapack.h"
#include "mex.h"
#include "qtcmn.h"

#define SQR(x)		((x)*(x))
#define SGN(x)		((x)<=0.0?-1.0:1.0)
#define ROUND(x)	(floor((x)+0.5))

/* C(mxn)=alpha * A(mxk) * B(kxn) + beta * C ---------------------------------*/
static void dmul(double beta, double *C, double alpha, const double *A,
				 const double *B, int m, int k, int n, char *tr)
{
	int lda=tr[0]=='T'?k:m,ldb=tr[1]=='T'?n:k;
	dgemm(tr,tr+1,&m,&n,&k,&alpha,(double *)A,&lda,(double *)B,&ldb,&beta,C,&m);
}
/* B(nxm)=A(nxn) \ B(nxm) ----------------------------------------------------*/
static int dsolv(const double *A, double *B, int n, int m, char *tr)
{
	int i,info,*ipiv=(int *)malloc(sizeof(int)*n);
	double *Aa=MAT(n,n);
	if (!ipiv||!Aa) return 0;
	COPY(A,n,n,Aa);
	dgetrf(&n,&n,Aa,&n,ipiv,&info); /* LU decomposition */
	if (!info) dgetrs(tr,&n,&m,Aa,&n,ipiv,B,&n,&info); /* solve linear-eq */
	free(ipiv); free(Aa);
	return !info;
}
/* LD factorization : Q=L'*diag(D)*L -----------------------------------------*/
static int LD(const double *Q, int n, double *L, double *D)
{
	int i,j,k;
	double a,*Qc=MAT(n,n);
	if (!Qc) mexErrMsgTxt("malloc error\n");
	
	COPY(Q,n,n,Qc);
	for (i=n-1;i>=0;i--) {
		if ((D[i]=Qc[i+n*i])<=0.0) {free(Qc); return 0;}
		a=sqrt(D[i]);
		for (j=0;j<=i;j++) L[i+n*j]=Qc[i+n*j]/a;
		for (j=0;j<=i-1;j++) {
			for (k=0;k<=j;k++) Qc[j+n*k]-=L[i+n*k]*L[i+n*j];
		}
		for (j=0;j<=i;j++) L[i+n*j]/=L[i+n*i];
	}
	free(Qc);
	return 1;
}
/* integer gauss transformation ----------------------------------------------*/
static void gauss(double *L, double *Z, int n, int i, int j)
{
	int k,mu;
	
	if ((mu=(int)ROUND(L[i+n*j]))==0) return;
	for (k=i;k<n;k++) L[k+n*j]-=(double)mu*L[k+n*i];
	for (k=0;k<n;k++) Z[k+n*j]-=(double)mu*Z[k+n*i];
}
/* permutations --------------------------------------------------------------*/
static void perm(double *L, double *D, int n, int j, double del, double *Z)
{
	int k;
	double eta,lam,a0,a1;
	
	eta=D[j]/del;
	lam=D[j+1]*L[j+1+n*j]/del;
	D[j]=eta*D[j+1];
	D[j+1]=del;
	for (k=0;k<=j-1;k++) {
		a0=L[j+n*k];
		a1=L[j+1+n*k];
		L[j+n*k]=-L[j+1+n*j]*a0+a1;
		L[j+1+n*k]=eta*a0+lam*a1;
	}
	L[j+1+n*j]=lam;
	for (k=j+2;k<n;k++) SWAP(L[k+n*j],L[k+n*(j+1)]);
	for (k=0;k<n;k++) SWAP(Z[k+n*j],Z[k+n*(j+1)]);
}
/* reduction (lambda decorrelation) : z=Z'*a, Qz=Z'*Q*Z=L'*diag(D)*L ---------*/
static void reduction(double *L, double *D, int n, double *Z)
{
	int i,j,k;
	double del;
	for (i=0;i<n;i++) Z[i+n*i]=1.0;
	j=n-2; k=n-2;
	while (j>=0) {
		/* integer gauss transformation */
		if (j<=k) for (i=j+1;i<n;i++) gauss(L,Z,n,i,j);
		
		del=D[j]+SQR(L[j+1+n*j])*D[j+1];
		if (del<D[j+1]) {
			perm(L,D,n,j,del,Z); /* permutations */
			k=j; j=n-2;
		}
		else j--;
	}
}
/* modified lambda search with shrinking ------------------------------------*/
static void search(const double *L, const double *D, const double *zs, int n,
				   int p, double *zn, double *fn)
{
	int i,j,k,nn=0,imax=0;
	double newdist,maxdist=1E99,y,*S,*dist,*zb,*z,*step;

	S=ZMAT(n,n); dist=VEC(n); zb=VEC(n); z=VEC(n); step=VEC(n);
	if (!S||!dist||!zb||!z||!step) mexErrMsgTxt("malloc error\n");
	
	k=n-1; dist[k]=0.0;
	zb[k]=zs[k];
	z[k]=ROUND(zb[k]); y=zb[k]-z[k]; step[k]=SGN(y);
	while (1) {
		newdist=dist[k]+SQR(y)/D[k];
		if (newdist<maxdist) {
			if (k!=0) {
				dist[--k]=newdist;
				for (i=0;i<=k;i++)
					S[k+n*i]=S[k+1+n*i]+(z[k+1]-zb[k+1])*L[k+1+n*i];
				zb[k]=zs[k]+S[k+n*k];
				z[k]=ROUND(zb[k]); y=zb[k]-z[k]; step[k]=SGN(y);
			}
			else {
				if (nn<p) {
					if (nn==0||newdist>fn[imax]) imax=nn;
					for (i=0;i<n;i++) zn[i+n*nn]=z[i];
					fn[nn++]=newdist;
				}
				else {
					if (newdist<fn[imax]) {
						for (i=0;i<n;i++) zn[i+n*imax]=z[i];
						fn[imax]=newdist;
						for (i=imax=0;i<p;i++) if (fn[imax]<fn[i]) imax=i;
					}
					maxdist=fn[imax];
				}
				z[0]+=step[0]; y=zb[0]-z[0]; step[0]=-step[0]-SGN(step[0]);
			}
		}
		else {
			if (k==n-1) break;
			else {
				k++;
				z[k]+=step[k]; y=zb[k]-z[k]; step[k]=-step[k]-SGN(step[k]);
			}
		}
	}
	for (i=0;i<p-1;i++) { /* sort by fn */
		for (j=i+1;j<p;j++) {
			if (fn[i]<fn[j]) continue;
			SWAP(fn[i],fn[j]);
			for (k=0;k<n;k++) SWAP(zn[k+n*i],zn[k+n*j]);
		}
	}
	free(S); free(dist); free(zb); free(z); free(step);
}
/* lambda integer least-square ----------------------------------------------*/
extern void Lambda(const double *a, const double *Q, int n, int nf, int flg,
				   double *af, double *f, double *Z)
{
	double *L=ZMAT(n,n),*D=VEC(n),*z=VEC(n);
	
	if (!L||!D||!z) mexErrMsgTxt("malloc error\n");
	
	if (LD(Q,n,L,D)) {
		if (flg) {
			reduction(L,D,n,Z);				/* reduction */
			dmul(0.0,z,1.0,Z,a,n,n,1,"TN");	/* z=Z'*a */
		}
		else COPY(a,n,1,z);
		search(L,D,z,n,nf,af,f);			/* modified lambda search */
		if (flg) {
			if (!dsolv(Z,af,n,nf,"T"))		/* a=Z'\z */
				mexPrintf("warning: lambda-back-trasformation error\n");
		}
	}
	else mexPrintf("warning: lambda-reduction error\n");
	
	free(L); free(D); free(z);
}
/* mex interface -------------------------------------------------------------*/
extern void mexFunction(int nargout, mxArray *argout[], int nargin,
					    const mxArray *argin[])
{
	int nf,na,flg=1;
	double *a,*Q,*af,*f,*Z;
	
	if (nargin<3) mexErrMsgTxt("argin error");
	if (nargout>3) mexErrMsgTxt("argout error"); 
	
	a=mxGetPr(argin[0]); na=mxGetM(argin[0]);
	Q=mxGetPr(argin[1]);
	nf=(int)*mxGetPr(argin[2]);
	if (nargin>=4) flg=(int)*mxGetPr(argin[3]);
	if (na!=mxGetM(argin[1])||na!=mxGetN(argin[1])) {
		mexErrMsgTxt("argin error");
	}
	if (na>0) {
		argout[0]=mxCreateDoubleMatrix(na,nf,mxREAL); af=mxGetPr(argout[0]);
		argout[1]=mxCreateDoubleMatrix(1,nf,mxREAL);  f=mxGetPr(argout[1]);
		if (nargout<3) Z=ZMAT(na,na);
		else {
			argout[2]=mxCreateDoubleMatrix(na,na,mxREAL); Z=mxGetPr(argout[2]);
		}
		Lambda(a,Q,na,nf,flg,af,f,Z);
		if (nargout<3) free(Z);
	}
	else {
		argout[0]=mxCreateDoubleMatrix(0,0,mxREAL);
		argout[1]=mxCreateDoubleMatrix(0,0,mxREAL);
		argout[2]=mxCreateDoubleMatrix(0,0,mxREAL);
	}
}
