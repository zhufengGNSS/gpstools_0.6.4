/*------------------------------------------------------------------------------
% [system] : GpsTools
% [module] : read dgrb gpv data
% [func]   : read dgrb gpv data
% [argin]  : file   = grib file path
% [argout] : data   = gpv data (struct array)
%                data(n) = dgrb data struct
%                  model : model id
%                           4 : JMA RSM
%                           5 : JMA MSM
%                  id    : parameter id
%                           1 : Pressure (hPa)
%                           4 : Temparature (C)
%                          13 : Humidity (%)
%                          23 : u-wind (m/s)
%                          24 : v-wind (m/s)
%                          42 : w-wind (hPa/s)
%                          49 : Rain (mm/s)
%                         102 : Geopotential Height (m)
%                         225 : Cloud
%                  press : pressure (hPa)
%                  time  : initial time of forecast (UTC)
%                          [year,month,day,hour,min,sec]
%                  ft    : forecast time (hr) [ft1,ft2,ft3,...]
%                  i     : x grid start/end index [i1,i2]
%                  j     : y grid start/end index [i1,i2]
%                  data  : dgrb data
%                      data(n,m,t) : grid (n-i1+1,m-j1+1) ft(t) data
% [note]   :
% [version]: $Revision: 16 $ $Date: 2008-12-12 15:49:30 +0900 (金, 12 12 2008) $
%            Copyright(c) 2004-2008 by T.Takasu, all rights reserved
% [history]: 04/05/13   0.1  new
%            08/11/30   0.2  suppress warning
%-----------------------------------------------------------------------------*/
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "mex.h"

#define NMAX		32768	/* max section counts in one file */
#define SECHLEN		44		/* section header length */

typedef struct {			/* sector data struct */
	int model;
	int id;
	int press;
	int time[6];
	int ft[3];
	int i[2],j[2];
	float *data;
} DGSEC;

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
/* find section header -------------------------------------------------------*/
static int findsec(FILE *fp, const char *head)
{
	int i,c,n=(int)strlen(head);
	char buff[16]="";
	while ((c=fgetc(fp))!=EOF) {
		for (i=1;i<n;i++) buff[i-1]=buff[i];
		buff[n-1]=(char)c;
		if (strcmp(buff,head)==0) return 1;
	}
	return 0;
}
/* read sector ---------------------------------------------------------------*/
static int readsec(FILE *fp, DGSEC *sec)
{
	static unsigned char bitmsk[]={0x80,0x40,0x20,0x10,0x08,0x04,0x02,0x01};
	unsigned char buff[64],*data,*p;
	int i,j,k,n,bit,slen,x,nx,ny,nt,nbit;
	double xx,ef,rf;
	if (fread(buff,1,SECHLEN,fp)!=SECHLEN) {
		mexPrintf("warning : dgrb format error : file trancated\n");
		return 0;
	}
	slen=getb(buff,2); 				/* section length */
	sec->model=buff[5];				/* model id */
	sec->id=buff[8];				/* type id */
	sec->press=getb(buff+10,2);		/* pressure */
	sec->time[0]=(buff[12]>=70?1900:2000)+buff[12];	/* time */
	for (i=1;i<=3;i++) sec->time[i]=buff[12+i];
	sec->time[4]=sec->time[5]=0;
	sec->ft[0]=buff[18];			/* forecast time start:end:interval */
	sec->ft[1]=buff[19];
	sec->ft[2]=buff[20]-200;
	sec->i[0]=getb(buff+24,2);		/* grid start index */
	sec->j[0]=getb(buff+26,2);
	sec->i[1]=getb(buff+28,2);		/* grid end index */
	sec->j[1]=getb(buff+30,2);
	nbit=getb(buff+32,2);			/* bit length */
	ef=pow(2.0,(double)(getsb(buff+34,2)));	/* e-factor */
	rf=(buff[36]&0x80?-1.0:1.0)*	/* r-factor */
	   (double)getb(buff+37,3)*pow(16.0,(double)((buff[36]&0x7F)-64))/16777216.0;
	if (!(data=(unsigned char *)malloc(slen-SECHLEN))) return 0;
	if (fread(data,1,slen-SECHLEN,fp)!=slen-SECHLEN) {
		mexPrintf("warning : dgrb format error : file trancated\n");
		free(data);
		return 0;
	}
	nx=sec->i[1]-sec->i[0]+1;
	ny=sec->j[1]-sec->j[0]+1;
	nt=(sec->ft[1]-sec->ft[0])/sec->ft[2]+1;
	if ((nx*ny*nt*nbit-1)/8+1!=slen-SECHLEN) {
		mexPrintf("warning : dgrb format error : section length\n");
		return 0;
	}
	if (!(sec->data=(float *)malloc(sizeof(float)*nx*ny*nt))) {
		free(data);
		return 0;
	}
	for (i=0,p=data,bit=0;i<nt;i++)
	for (j=0;j<nx;j++)
	for (k=0;k<ny;k++) {
		for (n=1,x=0;n<=nbit;n++,bit++) {
			if (bit>7) {bit=0; p++;}
			if (*p&bitmsk[bit]) x+=(1<<(nbit-n));
		}
		xx=(double)x; if (xx<-1E32||1E32<xx) xx=mxGetNaN();
		*(sec->data+i*nx*ny+j*ny+k)=(float)(rf+xx*ef);
	}
	free(data);
	return slen;
}
/* read dgrb file -----------------------------------------------------------*/
extern int readgrb(char *file, DGSEC *sec, int nmax)
{
	unsigned char buff[4];
	int nsec=0,len,slen;
	FILE *fp;
	
	if (!(fp=fopen(file,"rb"))) {
		mexPrintf("warning : dgrb file open error : %s\n",file);
		return 0;
	}
	while (findsec(fp,"DGRB")) {
		if (fread(buff,1,4,fp)!=4) break;
		for (len=getb(buff,2)-4;len>0;len-=slen,nsec++) {
			if (nsec>=nmax) {
				mexPrintf("warning : section counts overflow : %s\n",file);
				break;
			}
			if ((slen=readsec(fp,sec+nsec))==0) break;
		}
	}
	fclose(fp);
	return nsec;
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
extern void mexFunction(int nargout, mxArray *argout[], int nargin,
					    const mxArray *argin[])
{
	const char *fs[]={"model","id","press","time","ft","i","j","data"};
	static DGSEC sec[NMAX];
	static int id[NMAX],press[NMAX];
	char file[256];
	int i,j,n,m,t,nsec,ndata,nt,model,dims[3],ft[128];
	int imin=99999,imax=0,jmin=99999,jmax=0;
	float *x,*y;
	mxArray *p,*q;
	
	if (nargin<1||!mxIsChar(argin[0])) mexErrMsgTxt("argin error");
	if (mxGetString(argin[0],file,sizeof(file))!=0) mexErrMsgTxt("argin error");
	if (nargout>1) mexErrMsgTxt("argout error"); 
	
	/* read dgrb file */
	nsec=readgrb(file,sec,NMAX);
	
	if (nsec<=0) {
		argout[0]=mxCreateDoubleMatrix(0,0,mxREAL);
		return;
	}
	model=sec[0].model;
	for (i=ndata=0;i<nsec;i++) {
		for (j=0;j<6;j++)
			if (sec[i].time[j]!=sec[0].time[j])
				mexPrintf("warning : uncompatible time in %s\n",file);
		for (j=0;j<3;j++)
			if (sec[i].ft[j]!=sec[0].ft[j])
				mexPrintf("warning : uncompatible ft in %s\n",file);
		for (j=0;j<ndata;j++)
			if (id[j]==sec[i].id&&press[j]==sec[i].press) break;
		if (j>=ndata) {id[ndata]=sec[i].id; press[ndata++]=sec[i].press;}
		if (imin>sec[i].i[0]) imin=sec[i].i[0];
		if (imax<sec[i].i[1]) imax=sec[i].i[1];
		if (jmin>sec[i].j[0]) jmin=sec[i].j[0];
		if (jmax<sec[i].j[1]) jmax=sec[i].j[1];
	}
	for (nt=0;;nt++)
		if ((ft[nt]=sec[0].ft[0]+nt*sec[0].ft[2])>sec[0].ft[1]) break;
	
	p=mxCreateStructMatrix(ndata,1,8,fs);
	for (i=0;i<ndata;i++) {
		mxSetField(p,i,"model",CreateDouble(sec[0].model));
		mxSetField(p,i,"id",   CreateDouble(id[i]));
		mxSetField(p,i,"press",CreateDouble(press[i]));
		mxSetField(p,i,"time", CreateIntArray(sec[0].time,6));
		mxSetField(p,i,"ft",   CreateIntArray(ft,nt));
		dims[0]=imax-imin+1;
		dims[1]=jmax-jmin+1;
		dims[2]=nt;
		q=mxCreateNumericArray(3,dims,mxSINGLE_CLASS,mxREAL);
		x=(float *)mxGetPr(q);
		for (j=0;j<nsec;j++) {
			if (sec[j].id!=id[i]||sec[j].press!=press[i]) continue;
			for (t=0,y=sec[j].data;t<dims[2];t++)
			for (m=sec[j].j[0]-jmin;m<=sec[j].j[1]-jmin;m++)
			for (n=sec[j].i[0]-imin;n<=sec[j].i[1]-imin;n++)
				*(x+t*dims[0]*dims[1]+m*dims[0]+n)=*y++;
		}
		mxSetField(p,i,"data",q);
	}
	for (i=0;i<nsec;i++) free(sec[i].data);
	argout[0]=p;
}
