/*------------------------------------------------------------------------------
% [system] : GpsTools
% [module] : read grib gpv data
% [func]   : read grib gpv data
% [argin]  : file   = grib file path
% [argout] : data   = gpv data (struct array)
%                data(n).pprm = product parameter struct
%                  cid   : center id
%                  sid   : subcenter id
%                  pid   : process id
%                  gid   : grid id
%                          21-64 : international exchange grids (see GRIB specs)
%                          255 : non-defined grids
%                  param : parameter and units
%                            2 : mean sea level pressure(Pa)
%                            7 : geopotentiol height(m)
%                           11 : temperture(K)
%                           33 : wind u component(m/s)
%                           34 : wind v component(m/s)
%                           39 : vertical velocity(Pa/s)
%                           52 : relative humidity(%)
%                           81 : land ratio(1:land,0:sea,rate)
%                  level : layer level
%                          [  1,0] : surface
%                          [100,p] : pressure in p hPa
%                          [102,0] : mean sea level
%                  time  : initial time of forecast (UTC)
%                          [year,month,day,hour,min,sec]
%                  tunit : forecast time unit
%                            0 : min
%                            1 : hour
%                            2 : day
%                  p1    : period of time in tunit (0:analysis/initialized)
%                  p2    : period of time or time interval in tunit
%                data(n).gprm = grid parameter struct
%                  type  : data representation type
%                            0 : latitude/longitude grid (eq-cylindrical)
%                            1 : mercator projection
%                            3 : lambert conformal projection
%                            5 : polar stereographic projection
%                           13 : oblique lambert conformal projection
%                           50 : spherical harmonic coefficients
%                           90 : space view of orthgraphic grid
%                  nx(n) : no. of points along lon or y-axis
%                  ny    : no. of points along lat or x-axis
%                  rcflg : resolution and component flag
%                      rcflg(1) : direction increment given
%                      rcflg(2) : earth radius
%                            0 : sphere re = 6367.47km
%                            1 : oblate spheriod re=6378.16km f=1/297.0
%                      rcflg(3) : uv component of vector
%                            0 : relative to east/north
%                            1 : relative to defined grid x/y
%                  smode : scanning mode
%                            0 : points scan in + direction
%                            1 : points scan in - direction
%                  lat1,lon1 : lat/lon of first grid (deg or m)
%                  lat2,lon2 : lat/lon of last grid (deg) (type=0)
%                  dx,dy : lat/lon or x/y increment (deg or m)
%                  lov   : orientation of grid (y-parallel lon.) (deg) (type!=0)
%                  lati1,lati2 : lats which secant cone cuts the sphere (type!=0)
%                  latp,lonp : lat/lon of southern pole (type!=0)
%                data(n).data = grid data
%                  data(n).data(i,j) : grid (i,j) data (single precision)
% [note]   :
% [version]: $Revision: 4 $ $Date: 06/07/08 1:01 $
%            Copyright(c) 2004-2006 by T.Takasu, all rights reserved
% [history]: 04/05/05   0.1  new
%            04/05/13   0.2  single precision, search GRIB header
%-----------------------------------------------------------------------------*/
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "mex.h"

//#define DEBUG
#ifdef DEBUG
#define dmsg	mexPrintf
#else
#define dmsg	1?(void)0:(void)mexPrintf
#endif
#define emsg	mexPrintf
#define NDMAX	4096	/* max data counts in one file */
#define NXYMAX	512		/* max nx/ny */

typedef struct {		/* product parameters */
	int cid,sid,pid,gid;
	int param;
	int level[2];
	int time[6];
	int tunit;
	int p1,p2;
} PPRM;

typedef struct {		/* grid parameters */
	int type;
	int nxflg;
	int nx[NXYMAX];
	int ny;
	int rcflg[3];
	int smode[3];
	double lat1,lon1;
	double lat2,lon2;
	double dx,dy;
	double lov;
	double lati1,lati2;
	double latp,lonp;
} GPRM;

typedef struct {		/* grib data */
	PPRM pprm;			/* product parameters */
	GPRM gprm;			/* grid parameters */
	float *data;		/* grid data */
} GRIBD;

/* get unsigned bytes --------------------------------------------------------*/
extern int getb(unsigned char *p, int n)
{
	int x=0; for (;n>0;n--) x=(x<<8)+*p++; return x;
}
/* get signed bytes ----------------------------------------------------------*/
extern int getsb(unsigned char *p, int n)
{
	unsigned char buff[4];
	memcpy(buff,p,n); buff[0]&=0x7f;
	return (*p&0x80)?-getb(buff,n):getb(buff,n);
}
/* read unsigned bytes -------------------------------------------------------*/
static int readb(FILE *fp, int n, int *x)
{
	unsigned char buff[4];
	if (fread(buff,1,n,fp)!=(size_t)n) return 0;
	*x=getb(buff,n);
	return 1;
}
/* find grib header ----------------------------------------------------------*/
static int findhead(FILE *fp, const char *head)
{
	int i,c,n=strlen(head);
	char buff[16]="";
	while ((c=fgetc(fp))!=EOF) {
		for (i=1;i<n;i++) buff[i-1]=buff[i];
		buff[n-1]=(char)c;
		if (strcmp(buff,head)==0) return 1;
	}
	return 0;
}
/* read pds section ----------------------------------------------------------*/
static int readpds(FILE *fp, PPRM *pprm, int *flg, int *dscale)
{
	unsigned char buff[256];
	int i,length;
	if (!readb(fp,3,&length)||length<3+25||sizeof(buff)<length) {
		emsg("warning : grib pds length error : %d\n",length);
		return 0;
	}
	if (fread(buff,1,length-3,fp)!=length-3) {
		emsg("warning : grib pds trancated\n");
		return 0;
	}
	pprm->cid=buff[1];				/* center id */
	pprm->pid=buff[2];				/* process id */
	pprm->gid=buff[3];				/* grid id */
	flg[0]=(buff[4]&0x80)!=0;		/* gds flag */
	flg[1]=(buff[4]&0x40)!=0;		/* bms flag */
	pprm->param=buff[5];			/* parameter and units */
	pprm->level[0]=buff[6];			/* type */
	pprm->level[1]=getb(buff+7,2);	/* level */
	pprm->time[0]=buff[9]+(buff[9]>=70?1900:2000); /* initial time */
	for (i=1;i<5;i++) pprm->time[i]=buff[9+i];
	pprm->time[5]=0;
	pprm->tunit=buff[14];			/* forecast time unit */
	pprm->p1=buff[15];				/* period of time */
	pprm->p2=buff[16];				/* period of time or time interval */
	pprm->sid=buff[22];				/* subcenter id */
	*dscale=getsb(buff+23,2);		/* decimal factor (signed) */
	return 1;
}
/* read gds section ---------------------------------------------------------*/
static int readgds(FILE *fp, GPRM *gprm)
{
	unsigned char buff[256];
	int i,nv,pv,nx,dx,dy,length;
	if (!readb(fp,3,&length)||length<3+25||sizeof(buff)<length) {
		emsg("warning : grib gds length error : %d\n",length);
		return 0;
	}
	if (fread(buff,1,length-3,fp)!=length-3) {
		emsg("warning : grib gds trancated\n");
		return 0;
	}
	nv=buff[0];						/* vertical parameter number */
	pv=buff[1];						/* location of parameters list */
	gprm->type=buff[2];				/* data representation type */
	gprm->ny=getb(buff+5,2);		/* Nj or Ny */
	if (gprm->ny>NXYMAX) {
		emsg("warning : ny overflow : ny=%d\n",gprm->ny);
		return 0;
	}
	nx=getb(buff+3,2);				/* Ni or Nx */
	if (gprm->type==0&&nx==65535)
		for (i=0;i<gprm->ny;i++) gprm->nx[i]=getb(buff+pv-4+4*nv+i*2,2);
	else
		for (i=0;i<gprm->ny;i++) gprm->nx[i]=nx;
	gprm->lat1=(double)getsb(buff+7,3)*1E-3;  /* lat of 1st grid */ 
	gprm->lon1=(double)getsb(buff+10,3)*1E-3; /* lon of 1st grid */
	gprm->rcflg[0]=(buff[13]&0x80)!=0;
	gprm->rcflg[1]=(buff[13]&0x40)!=0;
	gprm->rcflg[2]=(buff[13]&0x08)!=0;
	gprm->smode[0]=(buff[24]&0x80)!=0;
	gprm->smode[1]=(buff[24]&0x40)!=0;
	gprm->smode[2]=(buff[24]&0x20)!=0;
	if (gprm->type==0) { /* lat/lon grid */
		gprm->lat2=(double)getsb(buff+14,3)*1E-3;
		gprm->lon2=(double)getsb(buff+17,3)*1E-3;
		dx=getb(buff+20,2);
		dy=getb(buff+22,2);
		gprm->dx=dx==65535||nx==65535?0.0:(double)dx*1E-3;
		gprm->dy=dy==65535?0.0:(double)dy*1E-3;
		gprm->lov=gprm->lati1=gprm->lati2=gprm->latp=gprm->lonp=0.0;
	}
	else {
		gprm->lov=(double)getsb(buff+14,3)*1E-3;
		gprm->dx=(double)getb(buff+17,3);
		gprm->dy=(double)getb(buff+20,3);
		gprm->lati1=(double)getsb(buff+25,3)*1E-3;
		gprm->lati2=(double)getsb(buff+28,3)*1E-3;
		gprm->latp=(double)getsb(buff+31,3)*1E-3;
		gprm->lonp=(double)getsb(buff+34,3)*1E-3;
		gprm->lat2=gprm->lon2=0.0;
	}
	return 1;
}
/* skip bms section ----------------------------------------------------------*/
static int readbms(FILE *fp)
{
	int i,length;
	if (!readb(fp,3,&length)) return 0;
	for (i=0;i<length-3;i++) if (fgetc(fp)==EOF) return 0;
	return 1;
}
/* read bds section ---------------------------------------------------------*/
static int readbds(FILE *fp, GRIBD *data, int dscale)
{
	static unsigned char bitmsk[]={0x80,0x40,0x20,0x10,0x8,0x4,0x2,0x1};
	unsigned char *buff;
	int i,j,x,nx,bit,off,length,gflg,ptype,dtype,aflg,ubits,bscale;
	int bchar,bmant,bsize;
	double bsign,ref,bfact,dfact;
	if (!readb(fp,3,&length)||length<3+8) {
		emsg("warning : grib bds length error : %d\n",length);
		return 0;
	}
	if (!(buff=(unsigned char *)malloc(length-3))) return 0;
	if (fread(buff,1,length-3,fp)!=length-3) {
	    free(buff);
		emsg("warning : grib bds trancated\n");
		return 0;
	}
	gflg =(buff[0]&0x80)!=0;	/* 0:gpv data */
	ptype=(buff[0]&0x40)!=0;	/* 0:simple pack,1:complex pack */
	dtype=(buff[0]&0x20)!=0;	/* 0:float,1:integer */
	aflg =(buff[0]&0x10)!=0;	/* 0:no add flag,1:add flag */
	ubits=buff[0]&0x0F;			/* unusebis */
	bscale=getsb(buff+1,2);		/* binary scale factor */
	bsign=(buff[3]&0x80)?-1.0:1.0; /* base sign */
	bchar=(buff[3]&0x7F)-64;	/* base charac */
	bmant=getb(buff+4,3);		/* base mantis */
	bsize=buff[7];				/* bit size */
	if (gflg!=0||ptype!=0||aflg!=0||bsize<=0||32<bsize) {
		emsg("warning : unsupported dbs param\n");
		free(buff);
		return 0;
	}
	ref=bsign*pow(2.0,-24.0)*(double)bmant*pow(16.0,(double)bchar);
	bfact=pow(2.0,(double)bscale);
	dfact=pow(10.0,(double)dscale);
	nx=((length-11)*8-ubits)/bsize;
	if (!(data->data=(float *)malloc(sizeof(float)*nx))) return 0;
	for (i=0,off=8,bit=0;i<nx;i++) {
		for (j=0,x=0;j<bsize;j++) {
			x=(x<<1)+(buff[off]&bitmsk[bit++]?1:0);
			if (bit>=8) {bit=0;off++;}
		}
		*(data->data+i)=(float)(ref+(double)x*bfact)/dfact;
	}
	free(buff);
	return 1;
}
/* read grib data -----------------------------------------------------------*/
extern int readgribdata(FILE *fp, GRIBD *data)
{
	int length,flg[2],dscale;
	char buff[8];
	if (!findhead(fp,"GRIB")) return 0;
	readb(fp,3,&length);
	fgetc(fp);
	if (!readpds(fp,&data->pprm,flg,&dscale)) {
		emsg("warning : grib pds error\n");
		return 0;
	}
	if (flg[0]&&!readgds(fp,&data->gprm)) {
		emsg("warning : grib gds error\n");
		return 0;
	}
	if (flg[1]&&!readbms(fp)) {
		emsg("warning : grib bms error\n");
		return 0;
	}
	if (!readbds(fp,data,dscale)) {
		emsg("warning : grib bds error\n");
		return 0;
	}
	fread(buff,1,4,fp); buff[4]='\0';
	if (strcmp(buff,"7777")!=0)
		emsg("warning : grib end section error\n");
	return 1;
}
/* read grib file -----------------------------------------------------------*/
extern int readgrib(char *file, GRIBD *data, int nmax)
{
	int n;
	FILE *fp;
	
	if (!(fp=fopen(file,"rb"))) {
		emsg("warning : grib file open error : %s\n",file);
		return 0;
	}
	for (n=0;n<nmax;n++) if (!readgribdata(fp,&data[n])) break;
	fclose(fp);
	return n;
}
/* mex interface -------------------------------------------------------------*/
static mxArray *CreateDouble(double x)
{
	mxArray *p=mxCreateDoubleMatrix(1,1,mxREAL);
	*mxGetPr(p)=x;
	return p;
}
static mxArray *CreateIntArray(const int *x, int n)
{
	int i;
	mxArray *p=mxCreateDoubleMatrix(1,n,mxREAL);
	for (i=0;i<n;i++) *(mxGetPr(p)+i)=(double)x[i];
	return p;
}
static mxArray *CreatePPrm(const PPRM *pprm)
{
	const char *fs[]={"cid","sid","pid","gid","param","level","time","tunit","p1","p2"};
	mxArray *p=mxCreateStructMatrix(1,1,10,fs);
	mxSetField(p,0,"cid",  CreateDouble(pprm->cid));
	mxSetField(p,0,"sid",  CreateDouble(pprm->sid));
	mxSetField(p,0,"pid",  CreateDouble(pprm->pid));
	mxSetField(p,0,"gid",  CreateDouble(pprm->gid));
	mxSetField(p,0,"param",CreateDouble(pprm->param));
	mxSetField(p,0,"level",CreateIntArray(pprm->level,2));
	mxSetField(p,0,"time", CreateIntArray(pprm->time,6));
	mxSetField(p,0,"tunit",CreateDouble(pprm->tunit));
	mxSetField(p,0,"p1",   CreateDouble(pprm->p1));
	mxSetField(p,0,"p2",   CreateDouble(pprm->p2));
	return p;
}
static mxArray *CreateGPrm(const GPRM *gprm)
{
	const char *fs[]={"type","nx","ny","rcflg","smode","lat1","lon1",
		"lat2","lon2","dx","dy","lov","lati1","lati2","latp","lonp"};
	mxArray *p;
	p=mxCreateStructMatrix(1,1,16,fs);
	mxSetField(p,0,"type", CreateDouble(gprm->type));
	mxSetField(p,0,"nx",   CreateDouble(gprm->nx[0]));
	mxSetField(p,0,"ny",   CreateDouble(gprm->ny));
	mxSetField(p,0,"rcflg",CreateIntArray(gprm->rcflg,3));
	mxSetField(p,0,"smode",CreateIntArray(gprm->smode,3));
	mxSetField(p,0,"lat1", CreateDouble(gprm->lat1));
	mxSetField(p,0,"lon1", CreateDouble(gprm->lon1));
	mxSetField(p,0,"lat2", CreateDouble(gprm->lat2));
	mxSetField(p,0,"lon2", CreateDouble(gprm->lon2));
	mxSetField(p,0,"dx",   CreateDouble(gprm->dx));
	mxSetField(p,0,"dy",   CreateDouble(gprm->dy));
	mxSetField(p,0,"lov",  CreateDouble(gprm->lov));
	mxSetField(p,0,"lati1",CreateDouble(gprm->lati1));
	mxSetField(p,0,"lati2",CreateDouble(gprm->lati2));
	mxSetField(p,0,"latp", CreateDouble(gprm->latp));
	mxSetField(p,0,"lonp", CreateDouble(gprm->lonp));
	return p;
}
static void InterpGrid(float *x, int nx, const float *y, int ny)
{
	int i;
	double n,dx=1.0/(double)(nx-1),dy=1.0/(double)(ny-1);
	for (i=0;i<nx;i++) {
		n=(double)i*dx/dy;
		x[i]=y[(int)n]*(floor(n)+1.0-n)+y[(int)n+1]*(n-floor(n));
	}
}
extern void mexFunction(int nargout, mxArray *argout[], int nargin,
					    const mxArray *argin[])
{
	static GRIBD data[NDMAX];
	const char *fs[]={"pprm","gprm","data"};
	char file[256];
	int i,j,k,n,nx,dims[2];
	float *x,*v,buff[NXYMAX];
	double dlon;
	mxArray *p,*q;
	
	if (nargin<1||!mxIsChar(argin[0])) mexErrMsgTxt("argin error");
	if (mxGetString(argin[0],file,sizeof(file))!=0) mexErrMsgTxt("argin error");
	if (nargout>1) mexErrMsgTxt("argout error"); 
	
	n=readgrib(file,data,NDMAX);
	p=mxCreateStructMatrix(n,1,3,fs);
	for (i=0;i<n;i++) {
		for (j=nx=0;j<data[i].gprm.ny;j++) {
			if (nx<data[i].gprm.nx[j]) nx=data[i].gprm.nx[j];
		}
		dims[0]=data[i].gprm.ny;
		dims[1]=nx;
		q=mxCreateNumericArray(2,dims,mxSINGLE_CLASS,mxREAL);
		for (j=0,x=(float *)mxGetPr(q),v=data[i].data;j<dims[0];j++) {
			if (fabs(data[i].gprm.dx)<1E-6) {
				InterpGrid(buff,nx,v,data[i].gprm.nx[j]);
				for (k=0;k<nx;k++) *(x+dims[0]*k+j)=buff[k];
			}
			else {
				for (k=0;k<nx;k++) *(x+dims[0]*k+j)=*(v+k);
			}
			v+=data[i].gprm.nx[j];
		}
		if (data[i].gprm.type==0&&fabs(data[i].gprm.dx)<1E-6) {
			data[i].gprm.nx[0]=nx;
			dlon=data[i].gprm.lon2-data[i].gprm.lon1;
			if (dlon<=-180) dlon+=360.0; else if (dlon>180) dlon-=360;
			data[i].gprm.dx=dlon/(nx-1);
		}
		mxSetField(p,i,"data",q);
		mxSetField(p,i,"pprm",CreatePPrm(&data[i].pprm));
		mxSetField(p,i,"gprm",CreateGPrm(&data[i].gprm));
		free(data[i].data);
	}
	argout[0]=p;
}
