/*-------------------------------------------------------------------------------
% [system] : GpsTools
% [module] : range model
% [func]   : ion-free phase pseudo-range model
% [argin]  : td    = date (mjd)
%            ts    = time (sec)
%            x,ix  = states, state indexes
%            iz    = obs. data index
%            stas  = satellite position/velocity (m,m/sec) (eci)
%            cdts  = satellite clock bias (m)
%            cdtr  = receiver clock bias (m)
%            zpdr  = zpdr(n) station n tropspheric zenith path delay (m)
%            posr  = station position (m) (ecef)
%            dpos  = station displacement (m) (ecef)
%            bcpr  = phase biases (m)
%            U     = eci to ecef convertion matrix
%            dx,dy,du = partial derivatives of U
%            rsun  = solar position (m) (ecef)
%            ants  = satellite antenna parameters [apc1;apc2;apv1;apv2] (38xn)
%                    apc1 = L1 phase center offset (m) [x;y;z] (sat fixed coordinate)
%                    apc2 = L2 phase center offset (m) [x;y;z] (sat fixed coordinate)
%                    apv1 = L1 phase center variation (m) (1x16)
%                    apv2 = L2 phase center variation (m) (1x16)
%                        apv?[i]=nadir(i) pcv (nadir=0:1:15deg)
%            antp  = receiver antenna parameters
%                        [apc1;apc2;ecc;apv1;apv2] (2783xm)
%                    apc1 = L1 phase center offset (m) [up;north;east]
%                    apc2 = L2 phase center offset (m) [up;north;east]
%                    ecc  = antenna deltas (m) [up;north;east]
%                    apv1 = L1 phase center variation (m) (73x19)
%                    apv2 = L2 phase center variation (m) (73x19)
%                        apv?[i+j*73]=az(i),el(j) pcv (az=0:5:360deg,el=0:5:90deg)
%            phs   = previous phase windup correction (rad)
%            elmin = min elevation angle (rad)
%            elmax = max elevation angle (rad)
%            metprm = meteological parameters
%            trop  = tropospheric model
%            mapf  = tropospheric mapping function
%            mpcc  = phase mutlpath parameters cnm
%            mpcs  = phase mutlpath parameters snm
%            f1,f2 = carrier frequencies L1,L2 (Hz)
%            obscorr = obs. corrections [apcs,apcr,rels,phw,mpr] (1:on)
%            clkref = reference clock flag (1:on)
%           (rmode) = receiver attitude model (0:fixed station,1:leo satellite)
%           (pmode) = receiver position model (0:ecef,1:eci)
%           (mapfc) = mapfc(:,n) tropospheric mapping function coefficients
% [argout] : g     = model (m)
%            dgds  = partial derivatives of model
%            azel  = satellite azimath/elevation angles (rad)
%            phn   = phase windup correction (rad)
%            index = valid obs data index
%           (gg)   = detailed model terms (for debug)
% [note]   :
% [version]: $Revision: 16 $ $Date: 2008-12-12 15:49:30 +0900 (金, 12 12 2008) $
%            Copyright(c) 2004-2006 by T.Takasu, all rights reserved
% [history]: 04/11/15  0.1  new
%            05/03/22  0.2  geocenter dx/dy/dz->dx/dy/dz/scale
%            05/05/31  0.2  add argout drds
%            05/06/26  0.3  add satellite antenna pcv correction
%            05/06/29  0.4  add argin trop,trgm
%            05/07/01  0.5  vectorized
%            05/08/11  0.6  add argin rmode
%            06/02/07  0.7  add argin bcpr,zpdr, argout gg for debug
%                           delete argin trgm
%            06/04/03  0.8  add argin pmode
%            06/05/07  0.9  add mapf_gmf as tropopheric mapping function
%            08/11/25  0.10 support no troposphere delay option (gt_0.6.4)
%                           support satellite L1/L2 antenna correction
%                           support receiver antenna azimuth angle dependancy
%                           support vmf1 for tropospheric mapping function
%                           add argin mapfc
%-----------------------------------------------------------------------------*/
#include <string.h>
#include "mex.h"
#include "qtcmn.h"

/*#define DEBUG*/

extern void EcefToGeod(const double *epos, double *gpos, double *E);
extern void EcsfToSatf(const double *state, double *E);
extern void SatRange(const double *state, const double *rpos, const double *U,
					 const double *dtr, int corrlightt, int direc, double *rs,
					 double *range, double *drds);
extern void SatAzEl(const double *spos, const double *rpos, double *azel);
extern void trop_saast(const double *t, const double *azel, const double *gpos,
					   const double *met_prm, double *tropd, double *tropw);
extern void mapf_nmf(double *t, const double *azel, const double *gpos,
					 double *mapfd, double *mapfw);
extern void SatApc(const double *rsat, const double *rrcv, const double *rsun,
				   const double *apc1, const double *apc2, const double *apv1,
				   const double *apv2, double *apcs);
extern void RcvApc(const double *azel, const double *apc1, const double *apc2,
				   const double *ecc, const double *apv1, const double *apv2,
				   double *apcr);
extern double RelCorr(const double *rsat, const double *vsat,
					  const double *rrcv, int opts);
extern double phwindup(const double *rsat, const double *rsun,
					   const double *posr, const double *U, const double *phwinp,
					   const double *rmode);
extern double RcvMpc(const double *azel, const double *mpcc, const double *mpcs,
					 int nmax);

extern void GMF(double *dmjd, double *dlat, double *dlon, double *dhgt,
				double *zd, double *gmfh, double *gmfw);
extern void VMF1_HT(double *ah, double *aw, double *dmjd, double *dlat,
				    double *dhgt, double *zd, double *vmf1h, double *vmf1w);

typedef struct {
	int ns,nr,nmpc;
	double *stas,*cdts,*cdtr,*zpdr,*posr,*dpos,*bcpr,*U,*dx,*dy,*du,*rsun;
	double *ants,*antp,*phs,*elmin,*elmax,*metr,*mpcc,*mpcs,*f1,*f2,cif[2];
	double *obscorr,*clkref,*rmode,*pmode,*mapfc;
	char trop[32],mapf[32];
} ModelP;

/* get state index (m,n:zero origin) -----------------------------------------*/
static int IxIndex(const mxArray *ix, const char *f, int m, int n, int *nx)
{
	mxArray *p;
	double *x;
	*nx=0;
	if ((p=mxGetField(ix,0,f))==NULL) return 0;
	if (mxIsCell(p)) {
		if ((int)mxGetM(p)<=m||(int)mxGetN(p)<=n||!(p=mxGetCell(p,m+(int)mxGetM(p)*n)))
			return 0;
	}
	if ((x=mxGetPr(p))==NULL) return 0;
	*nx=(int)mxGetN(p);
	return (int)x[0];
}
/* tropos. mapping function : COSZ -------------------------------------------*/
static void mapf_cosz(const double *t, const double *azel, const double *gpos,
					  double *mfd, double *mfw)
{
	*mfd=*mfw=1/sin(azel[1]);
}
/* eci state (t) to ecef position (t-dt) -------------------------------------*/
static void ecefpos(const double *state, const double *U, double dt,
				    double *pos)
{
	const double OMGE=7.2921151467E-5;
	double R[9],E[9],rsat[3],r=NORM(state);
	int i;
	for (i=0;i<3;i++) /* rs(t-dt) (eci) */
		rsat[i]=state[i]*(1.0+GME/2.0/r/r/r*dt*dt)-state[3+i]*dt;
	Rz(-OMGE*dt,R); MM(R,U,E); Mv(E,rsat,pos); /* eci->ecef */
}
/* range model ---------------------------------------------------------------*/
static int RangeModel(const double *td, const double *ts, const double *x,
				      int nx, const mxArray *ix, const double *iz,
				      const ModelP *p, double *g, double *dgds, double *azel,
				      double *phn, double *drds, double *gg)
{
	int i,j,s,r,n;
	double *p1,*p2,dt,stas[6],posr[6],cdtr,cdts,zpdr,dtr,rsat[3],spos[3],range;
	double gpos[3],Ut[9],rrcv[3],trd,trw,mfd=1.0,mfw=1.0,mfg[5],N;
	double ap[2],dp[3],dxdp[3],dydp[3],dudp[3],E[9],rs[3],rss[3];
	double t=td[0]+ts[0]/86400.0,zazel[]={0.0,M_PI/2.0},gx,gy;
	double c_apcs=0.0,c_apcr=0.0,c_rel=0.0,c_phw=0.0,c_mpc=0.0;
	double lat,lon,hgt,zd,*ah,*aw;
	
	dt=iz[0]-ts[0];			/* time offset */
	s=(int)iz[1]-1;			/* satellite index */
	r=(int)iz[2]-1;			/* station index */
	
	for (i=0;i<6;i++) {
		if (!mxIsNaN(p->stas[6*s+i])) continue;
#ifdef DEBUG
		mexPrintf("no valid sat orbit: t=%7,1f sat=%2d\n",*ts,s+1);
#endif
		return 0;
	}
	if (mxIsNaN(p->cdts[s])) {
#ifdef DEBUG
		mexPrintf("no valid sat clock: t=%7.1f sat=%2d\n",*ts,s+1);
#endif
		return 0;
	}
	if (mxIsNaN(p->cdtr[r])) {
#ifdef DEBUG
		mexPrintf("no valid rcv clock: t=%7.1f rcv=%2d\n",*ts,r+1);
#endif
		return 0;
	}
	for (i=0;i<6;i++) stas[i]=p->stas?p->stas[6*s+i]:0.0;
	for (i=0;i<6;i++) posr[i]=p->posr?p->posr[6*r+i]:0.0;
	cdts=p->cdts?p->cdts[s]:0.0;
	cdtr=p->cdtr?p->cdtr[r]-M_C*dt:0.0;
	zpdr=p->zpdr?p->zpdr[r]:mxGetNaN();

	if ((i=IxIndex(ix,"sato",0,s,&n))) {
		for (j=0;j<n;j++) stas[j]=x[i+j-1];
	}
	if ((i=IxIndex(ix,"rcvp",0,r,&n))) {
		for (j=0;j<n;j++) posr[j]=x[i+j-1];
	}
	if ((i=IxIndex(ix,"rcvc",0,r,&n))) cdtr=x[i-1];
	if ((i=IxIndex(ix,"rcvz",0,r,&n))) zpdr=x[i-1];
	dtr=p->cdtr[r]/M_C-dt;
	
	 /* receiver eci state to eci position */
	if ((int)p->pmode[0]==1) ecefpos(posr,p->U,dtr,posr);

	/* station displacement */
	for (i=0;i<3;i++) posr[i]+=p->dpos[3*r+i];
	
	/* geocenter offset */
	if ((i=IxIndex(ix,"eco",0,0,&n))) {
		for (j=0;j<n;j++) posr[j]=posr[j]*(1.0+x[i+2]*1E-9)+x[i+j-1];
	}
	/* satellite-station range and satellite tx position (eci) */
	SatRange(stas,posr,p->U,&dtr,1,1,rsat,&range,drds);
	
	/* satellite azimath/elevation angle */
	Mv(p->U,rsat,spos);
	
	if ((int)p->rmode[0]==0) /* fixed station */
	    SatAzEl(spos,posr,azel);
	else { /* leo satellite */
		EcsfToSatf(posr,E); /* eci->radial/alongtrk/crosstrk */
		for (i=0;i<3;i++) rs[i]=spos[i]-posr[i]; /* rcv->sat */
		Mv(E,rs,rss);
		azel[0]=atan2(-rss[2],rss[1]);
		azel[1]=asin(rss[0]/NORM(rss));
	}
	if (azel[1]<p->elmin[0]||p->elmax[0]<azel[1]) {
#ifdef DEBUG
		mexPrintf("elev out-of-range: t=%7.1f sat=%2d rcv=%2d el=%5.1f\n",
				  *ts,s+1,r+1,azel[1]*RAD2DEG);
#endif
		return 0;
	}
	/* tropospheric delay term */
	EcefToGeod(posr,gpos,NULL);
	if (strcmp(p->trop,"trop_saast")==0) {
		trop_saast(&t,zazel,gpos,p->metr+3*r,&trd,&trw);
	}
	else { /* no troposphere delay */
		trd=trw=0.0;
	}
	if (strcmp(p->mapf,"mapf_nmf")==0) { /* NMF */
		mapf_nmf(&t,azel,gpos,&mfd,&mfw);
	}
	else if (strcmp(p->mapf,"mapf_gmf")==0) { /* GMF */
		lat=gpos[0]*DEG2RAD; lon=gpos[1]*DEG2RAD;
		hgt=gpos[2]; zd=M_PI/2.0-azel[1];
		GMF(&t,&lat,&lon,&hgt,&zd,&mfd,&mfw);
	}
	else if (strcmp(p->mapf,"mapf_vmf1")==0&&p->mapfc) { /* VMF1 */
		ah=p->mapfc+2*r;
		aw=p->mapfc+2*r+1;
		lat=gpos[0]*DEG2RAD; lon=gpos[1]*DEG2RAD;
		hgt=gpos[2]; zd=M_PI/2.0-azel[1];
		VMF1_HT(ah,aw,&t,&lat,&hgt,&zd,&mfd,&mfw);
	}
	else { /* COSZ */
		mapf_cosz(&t,azel,gpos,&mfd,&mfw);
	}
	if (!mxIsNaN(zpdr)) trw=zpdr-trd;
	if ((i=IxIndex(ix,"rcvg",0,r,&n))) {
		gx=1.0/tan(azel[1])*cos(azel[0]);
		gy=1.0/tan(azel[1])*sin(azel[0]);
		mfg[0]=mfw*gx;
		mfg[1]=mfw*gy;
		mfw+=mfg[0]*x[i-1]+mfg[1]*x[i];
		if (n>=5) {
			mfg[2]=mfw*gx*gx;
			mfg[3]=mfw*gx*gy;
			mfg[4]=mfw*gy*gy;
			mfw+=mfg[2]*x[i+1]+mfg[3]*x[i+2]+mfg[4]*x[i+3];
		}
	}
	Tr(p->U,Ut); Mv(Ut,posr,rrcv);
	
	/* partial derivatives */
	for (i=0;i<nx;i++) dgds[i]=0.0;
	for (i=0;i<3;i++) dp[i]=-drds[i];
	if ((i=IxIndex(ix,"sato",0,s,&n))) for (j=0;j<6;j++) dgds[i+j-1]=drds[j];
	if ((i=IxIndex(ix,"satc",0,s,&n))) {dgds[i-1]=-1.0; cdts=x[i-1];}
	if ((i=IxIndex(ix,"rcvc",0,r,&n))) dgds[i-1]=1.0;
	else if ((int)p->clkref[r]) cdtr=0.0;
	if ((i=IxIndex(ix,"rcvz",0,r,&n))) dgds[i-1]=mfw;
	if ((i=IxIndex(ix,"rcvg",0,r,&n))) {
		for (j=0;j<2;j++) dgds[i+j-1]=mfg[j]*trw;
		if (n>=5) for (j=2;j<5;j++) dgds[i+j-1]=mfg[j]*trw;
	}
	if ((i=IxIndex(ix,"rcvp",0,r,&n))) {
		if ((int)p->pmode[0]==0) /* ecef postion */
			Mv(p->U,dp,dgds+i-1);
		else /* eci position */
			for (j=0;j<3;j++) dgds[i+j-1]=dp[j];
	}
	if ((i=IxIndex(ix,"erp",0,0,&n))) {
		Mv(p->dx,dp,dxdp); dgds[i-1]=DOT(dxdp,posr);
		Mv(p->dy,dp,dydp); dgds[i  ]=DOT(dydp,posr);
		Mv(p->du,dp,dudp); dgds[i+1]=DOT(dudp,posr);
	}
	if ((i=IxIndex(ix,"eco",0,0,&n))) {
		Mv(p->U,dp,dgds+i-1);
		for (j=0;j<3;j++) dgds[i+2]+=dgds[i-1+j]*1E-9*posr[j];
	}
	if ((i=IxIndex(ix,"arcn",s,r,&n))) {
		N=x[i-1]; dgds[i-1]=1.0;
	}
	else N=p->bcpr?p->bcpr[p->nr*s+r]:0.0;

	/* ion-free phase pseudo-rage model */
	*g=range+cdtr-cdts+mfd*trd+mfw*trw+N;
	
	/* correction terms */
	if ((int)p->obscorr[0]&&p->ants) {
		p1=p->ants+38*s;
		SatApc(rsat,rrcv,p->rsun,p1,p1+3,p1+6,p1+22,ap);
		c_apcs=p->cif[0]*ap[0]+p->cif[1]*ap[1];
	}
	if ((int)p->obscorr[1]&&p->antp) {
		p1=p->antp+2783*r;
		RcvApc(azel,p1,p1+3,p1+6,p1+9,p1+1396,ap);
		c_apcr=p->cif[0]*ap[0]+p->cif[1]*ap[1];
	}
	if ((int)p->obscorr[2]) {
		c_rel=RelCorr(rsat,stas+3,rrcv,(int)p->obscorr[2]>=2?1:0);
	}
	if ((int)p->obscorr[3]) {
		p1=p->phs+s+p->ns*r;
		p2=phn+s+p->ns*r;
		*p2=phwindup(rsat,p->rsun,posr,p->U,p1,p->rmode);
		c_phw=M_C*(p->cif[0]*(*p2)/p->f1[0]+p->cif[1]*(*p2)/p->f2[0])/2.0/M_PI;
	}
	if ((int)p->obscorr[4]&&p->mpcc&&p->mpcs) {
		p1=p->mpcc+(p->nmpc+1)*(p->nmpc+1)*r;
		p2=p->mpcs+(p->nmpc+1)*(p->nmpc+1)*r;
		c_mpc=RcvMpc(azel,p1,p2,p->nmpc);
	}
	*g+=c_apcs+c_apcr+c_rel+c_phw+c_mpc;
	
	if (gg) {
		gg[0]=posr[0]; gg[1]=posr[1]; gg[2]=posr[2]; gg[3]=dtr; gg[4]=range;
		gg[5]=cdtr; gg[6]=cdts; gg[7]=mfd*trd; gg[8]=mfw*trw; gg[9]=N;
		gg[10]=c_apcs; gg[11]=c_apcr; gg[12]=c_rel; gg[13]=c_phw; gg[14]=c_mpc;
	}
	if (mxIsNaN(*g)) {
#ifdef DEBUG
		mexPrintf("invalid model: t=%7.1f sat=%2d rcv=%2d (%d %d %d %d)\n",*ts,
                  s+1,r+1,mxIsNaN(c_apcs),mxIsNaN(c_apcr),mxIsNaN(c_rel),mxIsNaN(c_phw));
#endif
		return 0;
	}
	return 1;
}
/* vectorized range model ----------------------------------------------------*/
extern int rangemodel(const double *td, const double *ts, const double *x,
				      int nx, const mxArray *ix, const double *iz, int nz,
				      const ModelP *p, double *g, double *G, double *azel,
				      double *phn, double *drds, double *index, double *gg)
{
	int i,j,n;
	double izn[3];
	
	for (i=n=0;i<nz;i++) {
		for (j=0;j<3;j++) izn[j]=iz[i+nz*j];
		if (RangeModel(td,ts,x,nx,ix,izn,p,g+n,G+nx*n,azel+2*n,phn,drds+6*n,
					   gg?gg+15*i:NULL))
			index[n++]=(double)(i+1);
	}
	return n;
}
