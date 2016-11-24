/*-------------------------------------------------------------------------------
% [system] : GpsTools
% [module] : range model
% [argin]  : td    = date (mjd)
%            ts    = time (sec)
%            x,ix  = states, state indexes
%            iz    = iz(n,:) obs data n index [time,sat,rcv]
%            stas  = stas(n,:) satellite n position/velocity (m,m/sec) (eci)
%            cdts  = cdts(n) satellite n clock bias (m)
%            cdtr  = cdtr(n) station n receiver clock bias (m)
%            zpdr  = zpdr(n) station n tropspheric zenith path delay (m)
%            posr  = posr(:,n) station n position (m) (ecef)
%            dpos  = station displacement (m) (ecef)
%            bcpr  = phase biases (m)
%            U     = eci to ecef convertion matrix
%            dx,dy,du = partial derivatives of U
%            rsun  = solar position (m) (ecef)
%            ants  = ants(:,n) satellite n antenna offset [apso;apsv]
%                    apso = phase center offset (m) [x;y;z] (sat fixed coordinate)
%                    apsv = phase center variation (m) [na0;na1;na2;...;na15]
%                           (nadir angle 0-15deg)
%            antp  = antp(:,n) station n antenna parameters [apc1;apc2;ecc;apv1;apv2]
%                    apc1 = L1 phase center offset (m) [up;north;east]
%                    apc2 = L2 phase center offset (m) [up;north;east]
%                    ecc  = antenna deltas (m) [up;north;east]
%                    apv1 = L1 phase center variation (m) [el0;el5;el10;...;el90]
%                    apv2 = L2 phase center variation (m) [el0;el5;el10;...;el90]
%            phs   = phs(n,m) sat n - sta m previous phase windup correction (rad)
%            elmin = min elevation angle (rad)
%            elmax = max elevation angle (rad)
%            metr  = metr(:,n) station n meteological parameters [press;temp;humi]
%            trop  = tropospheric model
%            mapf  = tropospheric mapping function
%            mpcc  = mpcc(:,:,n) station n phase mutlpath parameters cnm
%            mpcs  = mpcs(:,:,n) station n phase mutlpath parameters snm
%            f1,f2 = carrier frequencies L1,L2 (Hz)
%            obscorr = obs. corrections [apcs,apcr,rels,phw,mpr] (1:on)
%            clkref = clkref(n) station n reference clock flag (1:on)
%           (rmode) = rcv attitude mode (0:fixed station,1:leo satellite)
%           (pmode) = rcv position mode (0:ecef,1:eci)
% [argout] : g     = g(n,1) obs index(n) model (m)
%            G     = G(n,:) obs index(n) partial derivatives of model by states
%            azel  = azel(n,:) obs index(n) satellite azimath/elevation angles (rad)
%            phs   = phs(n,m) sat n - sta m phase windup correction (rad)
%            drds  = drds(n,:) obs index(n) partial derivative of range by sat pos/vel
%            index = index(n) obs index corresponding to iz
%           (gg)   = detailed model terms (for debug)
%                    gg(n,1:3)= receiver postion (m)
%                    gg(n,4)  = sampling time offset (sec)
%                    gg(n,5)  = geometric range (m)
%                    gg(n,6)  = receiver clock (m)
%                    gg(n,7)  = satellite clock (m)
%                    gg(n,8)  = tropospheric delay (hydrostatic) (m)
%                    gg(n,9)  = tropospheric delay (wet) (m)
%                    gg(n,10) = phase bias (m)
%                    gg(n,11) = satellite antenna offset/pcv (m)
%                    gg(n,12) = receiver antenna offset/pcv (m)
%                    gg(n,13) = relativity correction (m)
%                    gg(n,14) = phase windup (m)
%                    gg(n,15) = phase multipath (m)
% [note]   :
% [version]: $Revision: 16 $ $Date: 06/07/20 14:26 $
%            Copyright(c) 2004-2006 by T.Takasu, all rights reserved
% [history]: 04/11/15  0.1  new
%            05/05/31  0.2  add argout drds
%            05/06/26  0.3  add satellite antenna pcv correction
%            05/06/29  0.4  add argin trop,trgm
%            05/07/01  0.5  vectorized
%            05/08/11  0.6  add argin rmode
%            06/02/07  0.7  add argin bcpr,zpdr, argout gg for debug
%                           delete argin trgm
%            06/04/03  0.8  add argin pmode
%-----------------------------------------------------------------------------*/
#include "mex.h"
#include "qtcmn.h"

typedef struct {
	int ns,nr,nmpc;
	double *stas,*cdts,*cdtr,*zpdr,*posr,*dpos,*bcpr,*U,*dx,*dy,*du,*rsun;
	double *ants,*antp,*elmin,*elmax,*metr,*mpcc,*mpcs,*f1,*f2,cif[2];
	double *obscorr,*clkref,*rmode,*pmode;
	char trop[32],mapf[32];
} ModelP;

extern int rangemodel(const double *td, const double *ts, const double *x,
				      int nx, const mxArray *ix, const double *iz, int nz,
				      const ModelP *p, double *g, double *G, double *azel,
				      double *phs, double *drds, double *index, double *gg);

/* get 3d size ---------------------------------------------------------------*/
static int Get3dSize(mxArray *arg, int *m)
{
	int i, n=mxGetNumberOfDimensions(arg),*ns=mxGetDimensions(arg);
	if (n<=1) {*m=0; return 0;}
	for (i=0;i<n;i++) {
		if (ns[i]<=0) continue;
		*m=n<=2?1:ns[2];
		return ns[0];
	}
	*m=0; return 0;
}
/* mex interface -------------------------------------------------------------*/
extern void mexFunction(int nargout, mxArray *argout[], int nargin,
					    const mxArray *argin[])
{
	const mxArray *ix;
	ModelP p;
	double *q,*td,*ts,*x,*iz,*phsp,*g,*G,*azel,*phs,*drds,*index,*gg=NULL,z0=0.0;
	int i,j,nx,nz,ns,nr,n,m,k;
	int argm[31]={1,1,1,0,1, 6,1,1,0,3, 3,0,3,3,3, 3,3,19,47,1, 1,1,3,0,0, 0,0,1,1,1, 1};
	int argn[31]={1,1,1,0,3, 1,1,1,0,1, 1,0,3,3,3, 3,1, 1, 1,1, 1,1,1,0,0, 0,0,1,1,5, 1};
	char msg[128];

	if (nargin<31) mexErrMsgTxt("argin error");
	nx=mxGetM(argin[2]);
	nz=mxGetM(argin[4]);
	p.ns=mxGetN(argin[5]);
	p.nr=mxGetM(argin[7]);
	argm[4]=nz;
	argn[5]=argm[6]=argn[17]=argm[19]=p.ns;
	argm[7]=argn[9]=argn[10]=argn[18]=argn[19]=argn[22]=p.nr;
	for (i=0;i<31;i++)
		if (mxGetM(argin[i])<argm[i]||mxGetN(argin[i])<argn[i]) {
			sprintf(msg,"argin(%d) error",i+1);
			mexErrMsgTxt(msg);
		}
	if (!mxIsStruct(argin[3])) mexErrMsgTxt("argin(4) error");
	if (!mxIsChar(argin[23])) mexErrMsgTxt("argin(24) error");
	if (!mxIsChar(argin[24])) mexErrMsgTxt("argin(25) error");
	n=Get3dSize(argin[25],&k); if (k!=0&&k!=p.nr) mexErrMsgTxt("argin(26) error");
	m=Get3dSize(argin[26],&k); if (k!=0&&k!=p.nr) mexErrMsgTxt("argin(27) error");
	if (n==m) p.nmpc=n-1; else mexErrMsgTxt("argin(26/27) error");
	if (nargout>7) mexErrMsgTxt("argout error"); 
	
	td    =mxGetPr(argin[0]);
	ts    =mxGetPr(argin[1]);
	x     =mxGetPr(argin[2]);
	ix    =argin[3];
	iz    =mxGetPr(argin[4]);
	p.stas=mxGetPr(argin[5]);
	p.cdts=mxGetPr(argin[6]);
	p.cdtr=mxGetPr(argin[7]);
	p.zpdr=mxGetPr(argin[8]);
	p.posr=mxGetPr(argin[9]);
	p.dpos=mxGetPr(argin[10]);
	p.bcpr=mxGetPr(argin[11]);
	p.U   =mxGetPr(argin[12]);
	p.dx  =mxGetPr(argin[13]);
	p.dy  =mxGetPr(argin[14]);
	p.du  =mxGetPr(argin[15]);
	p.rsun=mxGetPr(argin[16]);
	p.ants=mxGetPr(argin[17]);
	p.antp=mxGetPr(argin[18]);
	phs   =mxGetPr(argin[19]);
	p.elmin=mxGetPr(argin[20]);
	p.elmax=mxGetPr(argin[21]);
	p.metr=mxGetPr(argin[22]);
	mxGetString(argin[23],p.trop,sizeof(p.trop));
	mxGetString(argin[24],p.mapf,sizeof(p.mapf));
	p.mpcc=mxGetPr(argin[25]);
	p.mpcs=mxGetPr(argin[26]);
	p.f1  =mxGetPr(argin[27]);
	p.f2  =mxGetPr(argin[28]);
	p.cif[0]= p.f1[0]*p.f1[0]/(p.f1[0]*p.f1[0]-p.f2[0]*p.f2[0]);
	p.cif[1]=-p.f2[0]*p.f2[0]/(p.f1[0]*p.f1[0]-p.f2[0]*p.f2[0]);
	p.obscorr=mxGetPr(argin[29]);
	p.clkref=mxGetPr(argin[30]);
	p.rmode=nargin>31?mxGetPr(argin[31]):&z0;
	p.pmode=nargin>32?mxGetPr(argin[32]):&z0;
	
	g=VEC(nz); G=MAT(nx,nz); azel=MAT(2,nz); drds=MAT(6,nz); index=VEC(nz);
	if (!g||!G||!azel||!drds||!index) mexErrMsgTxt("malloc error");
	if (nargout>=7) {
		if (!(gg=MAT(15,nz))) mexErrMsgTxt("malloc error");
	}
	n=rangemodel(td,ts,x,nx,ix,iz,nz,&p,g,G,azel,phs,drds,index,gg);
	
	argout[0]=mxCreateDoubleMatrix(n, 1,mxREAL);
	argout[1]=mxCreateDoubleMatrix(n,nx,mxREAL);
	argout[2]=mxCreateDoubleMatrix(n, 2,mxREAL);
	argout[3]=mxCreateDoubleMatrix(p.ns,p.nr,mxREAL);
	argout[4]=mxCreateDoubleMatrix(n, 6,mxREAL);
	argout[5]=mxCreateDoubleMatrix(n, 1,mxREAL);
	
	COPY(g,n,1,mxGetPr(argout[0]));
	COPY(phs,p.ns,p.nr,mxGetPr(argout[3]));
	COPY(index,n,1,mxGetPr(argout[5]));
	for (i=0,q=mxGetPr(argout[1]);i<nx;i++) for (j=0;j<n;j++) *q++=*(G+nx*j+i);
	for (i=0,q=mxGetPr(argout[2]);i< 2;i++) for (j=0;j<n;j++) *q++=*(azel+2*j+i);
	for (i=0,q=mxGetPr(argout[4]);i< 6;i++) for (j=0;j<n;j++) *q++=*(drds+6*j+i);
	
	if (nargout>=7) {
		argout[6]=mxCreateDoubleMatrix(nz,15,mxREAL);
		for (i=0,q=mxGetPr(argout[6]);i<15;i++) for (j=0;j<nz;j++) *q++=*(gg+15*j+i);
	}
	free(g); free(G); free(azel); free(drds); free(index); free(gg);
}
