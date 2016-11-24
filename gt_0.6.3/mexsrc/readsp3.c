/*------------------------------------------------------------------------------
% [system] : GpsTools
% [system] : GpsTools
% [module] : read NGS sp3 ephemeris
% [func]   : read NGS sp3 format ephemeris
% [argin]  : file  = file name
% [argout] : epoch = start epoch [year,month,day,hour,min,sec]
%            time  = time vector (sec)
%            eph   = satellite ephemerides
%                    [x,y,z,cb(,vx,vy,vz,cd);...]
%                    x,y,z    : satellite position (ecef) (m)
%                    cb       : satellite clock bias (sec)
%                    vx,py,pz : satellite velocity (ecef) (m/sec)
%                    cd       : satellite clock drift (sec/sec)
%            sats  = satellite list
%                    ('GPSnn':GPS,'GLOnn':GLONASS,'LEOnn':LEO,'GALnn':GALILEO)
%           (accs) = accuracy of satellite orbits
%           (std)  = standard deviations for sp3c
%           (type) = file type
%                    ('':UNKNOWN,'G':GPS,'M':MIXED,'R':GLONASS,L:'LEO',E:'GALILEO)
%           (tsys) = time system
%                    ('':UNKNOWN,'GPS':GPS TIME,'UTC':UTC)
% [note]   :
% [version]: $Revision: 5 $ $Date: 06/07/08 1:01 $
%            Copyright(c) 2004-2006 by T.Takasu, all rights reserved
% [history]: 05/03/21  0.1  mfile->mex
%            05/05/21  0.2  support argouts accs,std,type,tsys and sp3c
%-----------------------------------------------------------------------------*/
#include <stdio.h>
#include <string.h>
#include "mex.h"
#include "qtcmn.h"

#define NSMAX		85		/* max satellite counts */

extern void CalToMjd(const double *dt, double *mjd, double *sec, int len_dt);

/* substring -----------------------------------------------------------------*/
static char *SubStr(const char *str, int pos, int len)
{
	static char buff[64];
	strncpy(buff,str+pos,len); buff[len]='\0';
	return buff;
}
/* string to number ----------------------------------------------------------*/
static double StrToNum(const char *str, int pos, int len)
{
	double value;
	if (strlen(str)<pos+len) return mxGetNaN();
	if (sscanf(SubStr(str,pos,len),"%lf",&value)==1) return value;
	return mxGetNaN();
}
/* sat code to sat name ------------------------------------------------------*/
static char *CodeToSat(char code)
{
	static char sat[2]="";
	if 	    (code=='R') return "GLO";
	else if (code=='L') return "LEO";
	else if (code=='E') return "GAL";
	else if (code=='G'||code==' ') return "GPS";
	sat[0]=code;
	return sat;
}
/* read sp3 header -----------------------------------------------------------*/
static int ReadHead(FILE *fp, double *epoch, char *type, char **sats,
				    double *accs, double *bfact, char *ftype, char *tsys)
{
	int i,j,no,ns=0,k=0,acc;
	char buff[256],ft[]="GMRLE";
	
	for (i=0;i<22;i++) {
		if (fgets(buff,sizeof(buff),fp)==NULL) break;
		
		if (i==0) {
			*type=buff[2];
			if (sscanf(SubStr(buff,3,28),"%lf %lf %lf %lf %lf %lf",
					   epoch,epoch+1,epoch+2,epoch+3,epoch+4,epoch+5)!=6)
				return 0;
		}
		else if (2<=i&&i<=6) {
			for (j=0;j<17;j++) {
				if ((no=(int)StrToNum(buff,10+3*j,2))<=0) continue;
				sprintf(sats[ns],"%s%02d",CodeToSat(buff[9+3*j]),no);
				accs[ns++]=mxGetNaN();
			}
		}
		else if (7<=i&&i<=11) {
			for (j=0;j<17;j++) {
				acc=(int)StrToNum(buff,10+3*j,2);
				if (acc>=1) accs[17*(i-7)+j]=pow(2.0,(double)acc)*1E-3;
				else accs[17*(i-7)+j]=mxGetNaN();
			}
		}
		else if (i==12) {
			strncpy(ftype,buff+3,1); ftype[1]='\0';
			strncpy(tsys,buff+9,3); tsys[3]='\0';
			for (j=0;j<5;j++) if (ftype[0]==ft[j]) break;
			if (j>=5) strcpy(ftype,"");
			if (strcmp(tsys,"GPS")!=0&&strcmp(tsys,"UTC")!=0) strcpy(tsys,"");
		}
		else if (i==14) {
			bfact[0]=StrToNum(buff,3,10);
			bfact[1]=StrToNum(buff,14,12);
		}
	}
	return ns;
}
/* read sp3 body -------------------------------------------------------------*/
static int ReadSp3Data(FILE *fp, double td, double ts, int nv, char **sats,
					   int ns, double *bfact, double **time, double **data,
					   double **std)
{
	double t[6],tdd,tss,value,fact[]={1E3,1E3,1E3,1E-6,0.1,0.1,0.1,1E-10};
	int i,j,k,nt=0,size=4096;
	char buff[256],sat[16];
	
	*time=(double *)malloc(sizeof(double)*size);
	*data=(double *)malloc(sizeof(double)*size*nv*ns);
	*std =(double *)malloc(sizeof(double)*size*nv*ns);
	if (*time==NULL||*data==NULL||*std==NULL) {
		free(*time); free(*data); free(*std);
		return 0;
	}
	for (i=0;i<size*nv*ns;i++) (*data)[i]=(*std)[i]=mxGetNaN();
	
	while (fgets(buff,sizeof(buff),fp)!=NULL) {
		if (!strcmp(SubStr(buff,0,3),"EOF")) break;
		
		if (sscanf(SubStr(buff,3,28),"%lf %lf %lf %lf %lf %lf",
				   t,t+1,t+2,t+3,t+4,t+5)!=6) {
			mexPrintf("warning : sp3 format error (epoch) : %s\n",buff);
			continue;
		}
		if (nt>=size) {
			size*=2;
			*time=(double *)realloc(*time,sizeof(double)*size);
			*data=(double *)realloc(*data,sizeof(double)*nv*ns*size);
			*std =(double *)realloc(*std, sizeof(double)*nv*ns*size);
			for (i=size*nv*ns/2;i<size*nv*ns;i++) (*data)[i]=(*std)[i]=mxGetNaN();
		}
		CalToMjd(t,&tdd,&tss,3);
		(*time)[nt]=(tdd-td)*86400.0+t[3]*3600.0+t[4]*60+t[5]-ts;
		for (i=0;i<ns*nv/4;i++) {
			if (fgets(buff,sizeof(buff),fp)==NULL) break;
			sprintf(sat,"%s%02.0f",CodeToSat(buff[1]),StrToNum(buff,2,2));
			for (j=0;j<ns;j++) if (strcmp(sats[j],sat)==0) break;
			if (j>=ns) {
				mexPrintf("warning : sp3 format error (sat) : %s\n",buff);
				continue;
			}
			for (k=0;k<4;k++) {
				value=StrToNum(buff,k*14+4,14);
				if (ABS(value-999999.999999)<1E-6) continue;
				if		(buff[0]=='P') (*data)[nt*ns*nv+j*nv+k]=value*fact[k];
				else if (buff[0]=='V') (*data)[nt*ns*nv+j*nv+k+4]=value*fact[k+4];
				else mexPrintf("warning : sp3 format error (type) : %s\n",buff);
				if (k<3) value=pow(bfact[0],StrToNum(buff,k*3+61,2))*1E-3;
				else     value=pow(bfact[1],StrToNum(buff,70,3))*1E-9;
				if		(buff[0]=='P') (*std)[nt*ns*nv+j*nv+k]=value;
				else if (buff[0]=='V') (*std)[nt*ns*nv+j*nv+k+4]=value*1E-4;
			}
		}
		nt++;
	}
	return nt;
}
/* read sp3 rinex clock data -----------------------------------------------------*/
extern int readsp3(const char *file, double *epoch, double **time, double **eph,
				   double **std, char **sats, double *accs, int *ns, int *nv,
				   char *ftype, char *tsys)
{
	FILE *fp;
	double td,ts,bfact[]={0.0,0.0};
	char type;
	int nt=0;
	if ((fp=fopen(file,"rt"))==NULL) return -1;
	*ns=ReadHead(fp,epoch,&type,sats,accs,bfact,ftype,tsys);
	*nv=type=='V'?8:4;
	CalToMjd(epoch,&td,&ts,6);
	if (*ns>0) nt=ReadSp3Data(fp,td,ts,*nv,sats,*ns,bfact,time,eph,std);
	fclose(fp);
	return nt;
}
/* mex interface -------------------------------------------------------------*/
extern void mexFunction(int nargout, mxArray *argout[], int nargin,
					    const mxArray *argin[])
{
	double epoch[6],*time,*eph,*std,*p,*q,*r, accs[NSMAX];
	char file[256],strs[NSMAX][16],*sats[NSMAX],ftype[4]="",tsys[8]="";
	int i,j,k,nt,ns,nv,dims[3];

	for (i=0;i<NSMAX;i++) sats[i]=strs[i];
	if (nargin<1||!mxIsChar(argin[0])||mxGetString(argin[0],file,sizeof(file))!=0)
		mexErrMsgTxt("argin error");
	if (nargout>8) mexErrMsgTxt("argout error"); 
	if ((nt=readsp3(file,epoch,&time,&eph,&std,sats,accs,&ns,&nv,ftype,tsys))<0) {
		for (i=0;i<nargout&&i<6;i++) argout[i]=mxCreateDoubleMatrix(0,0,mxREAL);
		for (i=6;i<nargout&&i<8;i++) argout[i]=mxCreateString("");
		return;
	}
	argout[0]=mxCreateDoubleMatrix(1,6,mxREAL); COPY(epoch,1,6,mxGetPr(argout[0]));
	argout[1]=mxCreateDoubleMatrix(nt,1,mxREAL); COPY(time,nt,1,mxGetPr(argout[1]));
	dims[0]=nt; dims[1]=nv; dims[2]=ns;
	argout[2]=mxCreateNumericArray(3,dims,mxDOUBLE_CLASS,mxREAL);
	p=mxGetPr(argout[2]);
	argout[3]=mxCreateCellArray(1,&ns);
	if (nargout>4) {
		argout[4]=mxCreateDoubleMatrix(ns,1,mxREAL);
		r=mxGetPr(argout[4]);
	}
	if (nargout>5) {
		argout[5]=mxCreateNumericArray(3,dims,mxDOUBLE_CLASS,mxREAL);
		q=mxGetPr(argout[5]);
	}
	if (nargout>6) argout[6]=mxCreateString(ftype);
	if (nargout>7) argout[7]=mxCreateString(tsys);
	for (i=0;i<nt;i++)
	for (j=0;j<ns;j++)
	for (k=0;k<nv;k++) {
		p[j*nt*nv+k*nt+i]=eph[i*ns*nv+j*nv+k];
		if (nargout>5) q[j*nt*nv+k*nt+i]=std[i*ns*nv+j*nv+k];
	}
	for (i=0;i<ns;i++) {
		mxSetCell(argout[3],i,mxCreateString(sats[i]));
		if (nargout>4) r[i]=accs[i];
	}
	free(time); free(eph); free(std);
}
