/*------------------------------------------------------------------------------
% [system] : GpsTools
% [module] : standard Kalman filter
% [func]   : calculate measurement update rule by standard Kalman filter
% [argin]  : x   = prefit states
%            P   = prefit covariences matrix
%            dz  = residuals
%            G   = jacobian matrix
%            R   = measurement noise corvariences matrix
%           (nsig) = outlier threshold (sigma,0:no exclusion)
% [argout] : xu  = postfit states
%            Pu  = postfit covarience matrix
%            vz  = valid residulas (1:valid,0:invalid/excluded)
%            stat= status (0:OK,0<>:error code)
% [note]   : use intel MKL CBLAS/LAPACK library
% [version]: $Revision: 12 $ $Date: 06/07/08 1:01 $
%            Copyright(c) 2004-2006 by T.Takasu, all rights reserved
% [history]: 04/07/22   0.1  new
%            05/03/09   0.2  use dpotrf(cholesky) instead of dgetrf(LU)
%                            add argout stat
%            06/07/04   0.3  return error code as stat, stat=0 means ok
%-----------------------------------------------------------------------------*/
#include "mkl_types.h"
#include "mkl_blas.h"
#include "mkl_lapack.h"
#include "mex.h"
#include "qtcmn.h"

/* C(mxn)=alpha * A(mxk) * B(kxn) + beta * C ---------------------------------*/
static void dmul(double beta, double *C, double alpha, const double *A,
				 const double *B, int m, int k, int n, char *tr)
{
	int lda=tr[0]=='T'?k:m,ldb=tr[1]=='T'?n:k;
	dgemm(tr,tr+1,&m,&n,&k,&alpha,(double *)A,&lda,(double *)B,&ldb,&beta,C,&m);
}
/* C(mxn)=alpha * A(mxn) * B(nxn) + beta * C (B:symmetric/upper tr) ----------*/
static void dmus(double beta, double *C, double alpha, const double *A,
				 const double *B, int m, int n)
{
	dsymm("R","U",&m,&n,&alpha,(double *)B,&n,(double *)A,&m,&beta,C,&m);
}
/* A(nxn)=A^-1 (A:symmetric positive-definite matrix/upper tr) ---------------*/
static int dins(double *A, int n)
{
	int info;
	if (n<=0) return -9999;
	dpotrf("U",&n,A,&n,&info); /* cholesky (upper tr) */
	if (!info) dpotri("U",&n,A,&n,&info); /* inversion (upper tr) */
	return info;
}
/* standard kalman filter ----------------------------------------------------*/
extern int FilterEkf(const double *x, const double *P, const double *dz,
				     const double *G, const double *R, double nsig, int nx,
				     int nz, double *xu, double *Pu, double *vz)
{
	double *F,*S,*Fj,*Sj,*dzj,*K,*p,*q;
	int i,j,k,*vv,nv,stat=0;

	COPY(x,nx,1,xu);
	COPY(P,nx,nx,Pu);
	if (nz<=0) return 0;
	
	F=MAT(nz,nx); S=MAT(nz,nz); vv=(int *)malloc(sizeof(int)*nz);
	if (!F||!S||!vv) mexErrMsgTxt("malloc error");
	
	COPY(R,nz,nz,S);
	dmus(0.0,F,1.0,G,P,nz,nx);							/* F=G*P */
	dmul(1.0,S,1.0,F,G,nz,nx,nz,"NT");					/* S=F*G'+R */
	
	for (i=nv=0;i<nz;i++) {
		if (nsig<=0.0||dz[i]*dz[i]<=S[nz*i+i]*nsig*nsig) vv[nv++]=i;
	}
	Fj=MAT(nx,nv); Sj=MAT(nv,nv); dzj=MAT(nv,1); K=MAT(nx,nv);
	if (!Fj||!Sj||!dzj||!K) mexErrMsgTxt("malloc error");
	
	for (i=0;i<nv;i++) {
		k=vv[i]; dzj[i]=dz[k];
		for (j=0,p=Fj+i*nx;j<nx;j++) *p++=F[k+j*nz];	/* Fj=F(vv,:)' */
		for (j=0,p=Sj+i*nv;j<=i;j++) *p++=S[k*nz+vv[j]]; /* Sj=S(vv,vv) */
	}
	if (!(stat=dins(Sj,nv))) { 							/* Sj=inv(Sj) */
		dmus(0.0,K,1.0,Fj,Sj,nx,nv);					/* K =Fj*Sj */
		dmul(1.0,xu,1.0,K,dzj,nx,nv,1,"NN");			/* xu=xu+K*dzj */
		dmul(1.0,Pu,-1.0,K,Fj,nx,nv,nx,"NT");			/* Pu=Pu-K*Fj' */
	}
	if (vz!=NULL) for (i=0;i<nv;i++) vz[vv[i]]=1.0;
	
	free(F); free(S); free(vv); free(Fj); free(Sj); free(dzj); free(K);
	return stat;
}
/* mex interface -------------------------------------------------------------*/
extern void mexFunction(int nargout, mxArray *argout[], int nargin,
					    const mxArray *argin[])
{
	int nz,nx,stat;
	double *x,*P,*G,*dz,*R,nsig=0.0,*xu,*Pu,*vz=NULL;

    if (nargin<5||mxGetM(argin[0])<1||mxGetN(argin[2])<1)
    	mexErrMsgTxt("argin error");
    if (nargout<2) mexErrMsgTxt("argout error"); 
    
    x =mxGetPr(argin[0]); nx=mxGetM(argin[0]);
    P =mxGetPr(argin[1]);
    dz=mxGetPr(argin[2]); nz=mxGetM(argin[2]);
    G =mxGetPr(argin[3]);
    R =mxGetPr(argin[4]);
    if (nargin>=6&&mxGetM(argin[5])>=1) nsig=*mxGetPr(argin[5]);
    
    if (nx!=mxGetM(argin[1])||nx!=mxGetN(argin[1])||nz!=mxGetM(argin[3])||
    	nx!=mxGetN(argin[3])||nz!=mxGetM(argin[4])||nz!=mxGetN(argin[4]))
    	mexErrMsgTxt("argin error");
	
	argout[0]=mxCreateDoubleMatrix(nx,1,mxREAL);  xu=mxGetPr(argout[0]);
	argout[1]=mxCreateDoubleMatrix(nx,nx,mxREAL); Pu=mxGetPr(argout[1]);
	if (nargout>2) {
		argout[2]=mxCreateDoubleMatrix(nz,1,mxREAL); vz=mxGetPr(argout[2]);
	}
	stat=FilterEkf(x,P,dz,G,R,nsig,nx,nz,xu,Pu,vz);
	if (nargout>3) {
		argout[3]=mxCreateDoubleMatrix(1,1,mxREAL);
		*mxGetPr(argout[3])=(double)stat;
	}
}
