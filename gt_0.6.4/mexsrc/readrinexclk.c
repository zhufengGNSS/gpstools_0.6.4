/*------------------------------------------------------------------------------
% [system] : GpsTools
% [module] : read rinex clock data
% [func]   : read rinex clock data file
% [argin]  : file   = file path
% [argout] : epoch  = first epoch [year,month,day,hour,min,sec]
%            time   = time vector relative to epoch
%            types  = clock data types
%            names  = satellite/station list
%            data   = clock data
%                     data(n,:) : index(n) clock data (sec)
%            index  = satellite/station index
%            sigs   = clock std. deviations
%                     sigs(n,:) : index(n) clock std. deviation (sec)
% [note]   :
% [version]: $Revision: 16 $ $Date: 2008-12-12 15:49:30 +0900 (金, 12 12 2008) $
%            Copyright(c) 2004-2008 by T.Takasu, all rights reserved
% [history]: 04/03/08   0.1  new
%            04/05/26   0.2  add argout sigs
%            06/02/16   0.3  suppress warning message
%            08/11/30   0.4  suppress warning
%-----------------------------------------------------------------------------*/
#include <stdio.h>
#include <string.h>
#include "mex.h"
#include "qtcmn.h"

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
	if ((int)strlen(str)<pos+len) return mxGetNaN();
	if (sscanf(SubStr(str,pos,len),"%lf",&value)==1) return value;
	return mxGetNaN();
}
/* unique string -------------------------------------------------------------*/
static int Unique(char *strs, int nstr, int size)
{
	int i,j;
	qsort(strs,nstr,size,(int(*)(const void *,const void *))strcmp);
	for (i=1,j=0;i<nstr;i++)
		if (strcmp(strs+j*size,strs+i*size)!=0)
			if (++j<i) strcpy(strs+j*size,strs+i*size);
	return j+1;
}
/* read rinex header ---------------------------------------------------------*/
static int ReadHead(FILE *fp, char *type, char **types)
{
	int i,ntype=0;
	char buff[256],*label;

	while (fgets(buff,sizeof(buff),fp)!=NULL) {
		label=buff+60;
		
		if (strstr(label,"RINEX VERSION / TYPE")!=NULL)
			*type=buff[20];
		else if (strstr(label,"# / TYPES OF DATA")!=NULL) {
			if (sscanf(buff,"%d",&ntype)!=1) continue;
			for (i=0;i<ntype;i++) {
				if (i!=0&&i%9==0) fgets(buff,sizeof(buff),fp);
				strcpy(types[i],SubStr(buff,10+(i%9)*6,2));
			}
		}
		else if (strstr(label,"END OF HEADER")!=NULL) break;
	}
	return ntype;
}
/* read rinex clock body -----------------------------------------------------*/
static int ReadData(FILE *fp, int ntype, double *epoch, double **time,
					char **names, double **data, double **sigs, double **index,
					int *nname)
{
	double t,td0=0.0,ts0=0.0,td,ts,ep[6];
	int i,j,nz=0,size=4096,prn;
	char buff[256],*strs;
	
	*time =(double *)malloc(size*sizeof(double));
	*index=(double *)malloc(size*sizeof(double));
	*data =(double *)malloc(size*sizeof(double));
	*sigs =(double *)malloc(size*sizeof(double));
	strs=(char *)malloc(size*8*sizeof(char));
	if (*time==NULL||*index==NULL||*data==NULL||*sigs==NULL||strs==NULL) return 0;
	
	while (fgets(buff,sizeof(buff),fp)!=NULL) {
		if (sscanf(SubStr(buff,8,26),"%lf %lf %lf %lf %lf %lf",&ep[0],&ep[1],
				   &ep[2],&ep[3],&ep[4],&ep[5])!=6) continue;
		CalToMjd(ep,&td,&ts,6);
		if (td0==0.0) {
			for (i=0;i<6;i++) epoch[i]=ep[i];
			td0=td; ts0=ts;
		}
		t=(td-td0)*86400.0+ts-ts0;
		if (nz>=size) {
			size*=2;
			*time =(double *)realloc(*time, size*sizeof(double));
			*index=(double *)realloc(*index,size*sizeof(double));
			*data =(double *)realloc(*data, size*sizeof(double));
			*sigs =(double *)realloc(*sigs, size*sizeof(double));
			strs=(char *)realloc(strs,size*8*sizeof(char));
			if (*time==NULL||*index==NULL||*data==NULL||*sigs==NULL||strs==NULL)
				return 0;
		}
		(*time)[nz]=t;
		strcpy(strs+8*nz,SubStr(buff,3,4));
		(*data)[nz]=StrToNum(buff,40,20);
		(*sigs)[nz]=StrToNum(buff,60,20);
		nz++;
	}
	*names=(char *)malloc(nz*8*sizeof(char));
	memcpy(*names,strs,nz*8*sizeof(char));
	*nname=Unique(*names,nz,8);
	for (i=0;i<nz;i++) {
		for (j=0;j<*nname;j++) {
			if (strcmp(strs+8*i,*names+8*j)!=0) continue;
			(*index)[i]=j+1;
			break;
		}
	}
	for (i=0;i<*nname;i++)
		if (sscanf(*names+8*i,"G%d",&prn)==1) sprintf(*names+8*i,"GPS%02d",prn);
	free(strs);
	return nz;
}
/* read rinex clock data -----------------------------------------------------*/
extern int ReadRinexClk(const char *file, double *epoch, double **time,
					    char **types, char **names, double **data, double **sigs,
					    double **index, int *ntype, int *nname)
{
	FILE *fp;
	char type;
	int nz;
	if ((fp=fopen(file,"rt"))==NULL) return -1;
	*ntype=ReadHead(fp,&type,types);
	if (*ntype<=0||type!='C') {
		fclose(fp);
		return -1;
	}
	nz=ReadData(fp,*ntype,epoch,time,names,data,sigs,index,nname);
	fclose(fp);
	return nz;
}
/* mex interface -------------------------------------------------------------*/
extern void mexFunction(int nargout, mxArray *argout[], int nargin,
					    const mxArray *argin[])
{
	double epoch[6],*time,*data,*sigs,*index;
	char file[256],strs[16][16],*types[16],*names;
	int i,nz,ntype,nname;

	for (i=0;i<16;i++) types[i]=strs[i];
	if (nargin<1||!mxIsChar(argin[0])) mexErrMsgTxt("argin error");
	if (mxGetString(argin[0],file,sizeof(file))!=0) mexErrMsgTxt("argin error");
	if (nargout>7) mexErrMsgTxt("argout error"); 
	
    nz=ReadRinexClk(file,epoch,&time,types,&names,&data,&sigs,&index,&ntype,
    				&nname);
	if (nz<0) {
		argout[0]=mxCreateDoubleMatrix(0,0,mxREAL);
		argout[1]=mxCreateDoubleMatrix(0,0,mxREAL);
		ntype=0;
		argout[2]=mxCreateCellArray(0,&ntype);
		argout[3]=mxCreateCellArray(0,&ntype);
		argout[4]=mxCreateDoubleMatrix(0,0,mxREAL);
		argout[5]=mxCreateDoubleMatrix(0,0,mxREAL);
		argout[6]=mxCreateDoubleMatrix(0,0,mxREAL);
		return;
	}
	argout[0]=mxCreateDoubleMatrix(1,6,mxREAL); COPY(epoch,1,6,mxGetPr(argout[0]));
	argout[1]=mxCreateDoubleMatrix(nz,1,mxREAL); COPY(time,nz,1,mxGetPr(argout[1]));
	argout[2]=mxCreateCellArray(1,&ntype);
	for (i=0;i<ntype;i++) mxSetCell(argout[2],i,mxCreateString(types[i]));
	argout[3]=mxCreateCellArray(1,&nname);
	for (i=0;i<nname;i++) mxSetCell(argout[3],i,mxCreateString(names+8*i));
	argout[4]=mxCreateDoubleMatrix(nz,1,mxREAL); COPY(data, nz,1,mxGetPr(argout[4]));
	argout[5]=mxCreateDoubleMatrix(nz,1,mxREAL); COPY(index,nz,1,mxGetPr(argout[5]));
	argout[6]=mxCreateDoubleMatrix(nz,1,mxREAL); COPY(sigs, nz,1,mxGetPr(argout[6]));
	free(time); free(names); free(data); free(sigs); free(index);
}
