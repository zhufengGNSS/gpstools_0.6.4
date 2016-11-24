/*------------------------------------------------------------------------------
% [system] : GpsTools
% [module] : square root convergence filter
% [func]   : calculate measurement update rule by square root convergence filter
% [argin]  : x   = prefit states
%            P   = prefit covariences matrix
%            dz  = residuals
%            G   = jacobian matrix
%            R   = measurement noise corvariences matrix
%           (nsig) = outlier threshold (sigma,0:no exclusion)
% [argout] : xu  = postfit states
%            Pu  = postfit covarience matrix
%            vz  = valid residulas (1:valid,0:invalid/excluded)
%            stat= status (1:OK,0:error)
% [note]   : use intel MKL CBLAS/LAPACK library
% [version]: $Revision: 12 $ $Date: 2008-11-25 10:02:15 +0900 (火, 25 11 2008) $
%            Copyright(c) 2004-2006 by T.Takasu, all rights reserved
% [history]: 04/07/22   0.1  new
%            05/07/22   0.2  add argout stat
%-----------------------------------------------------------------------------*/
#include "mkl_types.h"
#include "mkl_cblas.h"
#include "mkl_lapack.h"
#include "mex.h"
#include "qtcmn.h"

#define _NN			0
#define _TN			1
#define _NT			2
#define _TT			3

/* matrix product by cblas ---------------------------------------------------*/
static void BlasMul(double alpha, const double *A, const double *B, int tr,
				    int m, int k, int n, double beta, double *C)
{
	CBLAS_TRANSPOSE transa,transb;
	int lda,ldb,ldc=m;
	transa=(tr&1)?CblasTrans:CblasNoTrans;
	transb=(tr&2)?CblasTrans:CblasNoTrans;
	lda=(tr&1)?k:m;
	ldb=(tr&2)?n:k;
	cblas_dgemm(CblasColMajor,transa,transb,m,n,k,alpha,A,lda,B,ldb,beta,C,ldc);
}
/* matrix inversion by lapack ------------------------------------------------*/
static int BlasInv(double *A, int n)
{
	int *ipiv,lwork=64*n,info=1;
	double *work=VEC(lwork);
	
	ipiv=(int *)malloc(sizeof(int)*n);
	
	if (work&&ipiv) {
		dgetrf(&n,&n,A,&n,ipiv,&info);
		if (!info) dgetri(&n,A,&n,ipiv,work,&lwork,&info);
	}
	free(work); free(ipiv);
	return !info;
}
/* cholesky decomposition ----------------------------------------------------*/
static int Chol(double *A, int n)
{
	char uplo='L';
	int i,j,info=1;
	dpotrf(&uplo,&n,A,&n,&info); /* L=chol(A)'; */
	if (!info)
		for (i=1;i<n;i++) for (j=0;j<i;j++) A[i*n+j]=0.0; /* clear upper tr */
	return !info;
}
/* square root convergence filter --------------------------------------------*/
extern int FilterSrcf(const double *x, const double *P, const double *dz,
					   const double *G, const double *R, double nsig, int nx,
					   int nz, double *xu, double *Pu, double *vz)
{
	double *C,*V,*S,*Vj,*Sj,*dzj,*K,*VS,*W,*Cu;
	int i,j,*vv,nvz,stat=0;

	COPY(x,nx,1,xu);
	COPY(P,nx,nx,Pu);
	if (nz<=0) return 1;
	
	C=MAT(nx,nx); V=MAT(nx,nz); S=MAT(nz,nz); vv=(int *)malloc(sizeof(int)*nz);
	if (!C||!V||!S||!vv) mexErrMsgTxt("malloc error");
	
	COPY(P,nx,nx,C);
	if (!Chol(C,nx)) {									/* C=chol(P)' */
		FreeMat(C); FreeMat(V); FreeMat(S); free(vv);
		COPY(P,nx,nx,Pu);
		return 0;
	}
	BlasMul(1.0,C,G,_TT,nx,nx,nz,0.0,V);				/* V=C'*G' */
	COPY(R,nz,nz,S);									/* S=V'*V+R */
	BlasMul(1.0,V,V,_TN,nz,nx,nz,1.0,S);
	
	for (i=nvz=0;i<nz;i++) {
		if (nsig<=0.0||dz[i]*dz[i]<=S[nz*i+i]*nsig*nsig) vv[nvz++]=i;
	}
	Vj=MAT(nx,nvz); Sj=MAT(nvz,nvz); dzj=MAT(nvz,1); K=MAT(nx,nvz);
	VS=MAT(nx,nvz); W=ZMAT(nx,nx); Cu=MAT(nx,nx);
	if (!Vj||!Sj||!dzj||!K||!VS||!W||!Cu) mexErrMsgTxt("malloc error");
	
	for (i=0;i<nvz;i++) {
		dzj[i]=dz[vv[i]];
		for (j=0;j<nx;j++) {
			Vj[i*nx+j]=V[vv[i]*nx+j];
		}
		for (j=0;j<nvz;j++) Sj[i*nvz+j]=S[vv[i]*nz+vv[j]];
	}
	if ((stat=BlasInv(Sj,nvz))) { 						/* Sj=inv(Sj) */
		BlasMul(1.0,Vj,Sj, _NN,nx,nvz,nvz,0.0,VS);		/* VS=Vj*Sj */
		BlasMul(1.0,C, VS, _NN,nx,nx, nvz,0.0,K );		/* K =C*VS */
		BlasMul(1.0,K, dzj,_NN,nx,nvz,1,  1.0,xu);		/* xu=x+K*dzj */
		for (i=0;i<nx;i++) W[i*nx+i]=1.0;				/* W =I-VS*Vj' */
		BlasMul(-1.0,VS,Vj,_NT,nx,nvz,nx,1.0,W);
		if ((stat=Chol(W,nx))) {						/* Cu=C*chol(W)' */
			BlasMul(1.0,C, W, _NN,nx,nx,nx,0.0,Cu);
			BlasMul(1.0,Cu,Cu,_NT,nx,nx,nx,0.0,Pu);		/* Pu=Cu*Cu' */
		}
	}
	if (vz!=NULL) for (i=0;i<nvz;i++) vz[vv[i]]=1.0;
	
	FreeMat(C); FreeMat(V); FreeMat(S); free(vv); FreeMat(Vj); FreeMat(Sj);
	FreeMat(dzj); FreeMat(K); FreeMat(VS); FreeMat(W); FreeMat(Cu);
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
	stat=FilterSrcf(x,P,dz,G,R,nsig,nx,nz,xu,Pu,vz);
	if (nargout>3) {
		argout[3]=mxCreateDoubleMatrix(1,1,mxREAL);
		*mxGetPr(argout[3])=(double)stat;
	}
}
