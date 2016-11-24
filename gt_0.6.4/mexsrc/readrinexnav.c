/*------------------------------------------------------------------------------
% [system] : GpsTools
% [module] : read rinex navigation message
% [func]   : read rinex navigation message file
% [argin]  : file    = file path
% [argout] : sats    = satellites list
%            rcv     = receiving station
%            data    = navigation messages
%                      data(n,:) : index(n) navigation message
%            index   = satellite index
%            ionprm  = ionospheric parameters
%                      ionprm(1,:) : ion alpha
%                      ionprm(2,:) : ion beta
%            dutc    = delta utc parameters [A0,A1,T,W]
%            comment = rinex header separated by ';'
% [note]   :
% [version]: $Revision: 16 $ $Date: 2008-12-12 15:49:30 +0900 (金, 12 12 2008) $
%            Copyright(c) 2004-2008 by T.Takasu, all rights reserved
% [history]: 04/11/17  0.1  new
%            08/11/30  0.2  suppress warning
%-----------------------------------------------------------------------------*/
#include <stdio.h>
#include <string.h>
#include "mex.h"
#include "qtcmn.h"

#define SIZE_NAV		37			/* record size */

/* substring -----------------------------------------------------------------*/
static char *SubStr(const char *str, int pos, int len)
{
	static char buff[64];
	if ((int)strlen(str)<pos+len) return "";
	strncpy(buff,str+pos,len); buff[len]='\0';
	return buff;
}
/* string to number ----------------------------------------------------------*/
static double StrToNum(char *str, int pos, int len)
{
	double value;
	char *p;
	for (p=str;(p=strchr(p,'D'))!=NULL;p++) *p='E';
	if (sscanf(SubStr(str,pos,len),"%lf",&value)==1) return value;
	return mxGetNaN();
}
/* read rinex header ---------------------------------------------------------*/
static void ReadHead(FILE *fp, char *type, double *ionprm, double *dutc,
				     char *comment)
{
	char buff[256],*label,*p=comment;
	int i;

	while (fgets(buff,sizeof(buff),fp)!=NULL) {
		strcpy(p,buff); p+=strlen(p)-1; strcpy(p++,";");
		label=buff+60;
		
		if (strstr(label,"RINEX VERSION / TYPE")!=NULL) {
			*type=buff[20];
		}
		else if (strstr(label,"ION ALPHA")!=NULL) {
			for (i=0;i<4;i++) ionprm[i*2]=StrToNum(buff,3+i*12,11);
		}
		else if (strstr(label,"ION BETA")!=NULL) {
			for (i=0;i<4;i++) ionprm[i*2+1]=StrToNum(buff,3+i*12,11);
		}
		else if (strstr(label,"DELTA-UTC:")!=NULL) {
			for (i=0;i<2;i++) dutc[i]=StrToNum(buff,3+i*19,19);
			for (i=0;i<2;i++) dutc[i+2]=StrToNum(buff,41+i*9,9);
		}
		else if (strstr(label,"END OF HEADER")!=NULL) break;
	}
}
/* read rinex navigation message body ----------------------------------------*/
static int ReadNav(FILE *fp, int nline, double **data, double **index)
{
	int i,nnav=0,line=0,size=1024;
	char buff[256];
	double *p;
	
	*index=(double *)malloc(size*sizeof(double));
	*data=(double *)malloc(size*SIZE_NAV*sizeof(double)); p=*data;
	if (*index==NULL||*data==NULL) return 0;
	
	while (fgets(buff,sizeof(buff),fp)!=NULL) {
		line=line%nline+1;
		if (line==1) {
			*((*index)+nnav)=0.0;
			sscanf(SubStr(buff,0,2),"%lf", (*index)+nnav);
			for (i=0;i<SIZE_NAV;i++) *(p+i)=0.0;
			sscanf(SubStr(buff,2,30),"%lf %lf %lf %lf %lf %lf",p,p+1,p+2,p+3,p+4,p+5);
			for (i=1;i<4;i++) *(p+5+i)=StrToNum(buff,i*19+3,19);
		}
		else {
			for (i=0;i<4;i++) *(p+line*4+i+1)=StrToNum(buff,i*19+3,19);
		}
		if (line==nline) {
			if (nnav>=size) {
				*index=(double *)realloc(*index,size*sizeof(double));
				*data=(double *)realloc(*data,size*SIZE_NAV*sizeof(double));
			}
			nnav++;
			p=*data+nnav*SIZE_NAV;
		}
	}
	return nnav;
}
/* read rinex navigation message file ----------------------------------------*/
extern int ReadRinexNav(const char *file, double **data, double **index,
						double *ionprm, double *dutc, char *comment, char *type)
{
	FILE *fp;
	int nnav=0;
	
	if ((fp=fopen(file,"rt"))==NULL) return -1;
	ReadHead(fp,type,ionprm,dutc,comment);
	if		(*type=='N') nnav=ReadNav(fp,8,data,index);
	else if (*type=='G') nnav=ReadNav(fp,4,data,index);
	else {
		fclose(fp);
		return 0;
	}
	fclose(fp);
	return nnav;
}
/* mex interface -------------------------------------------------------------*/
extern void mexFunction(int nargout, mxArray *argout[], int nargin,
					    const mxArray *argin[])
{
	double *data,*index,ionprm[8],dutc[4],*data_b;
	char *p,type,file[256],comment[4096],sat[16],rcv[16];
	int i,j,k,nnav,nsat=31;
	
	if (nargin<1||!mxIsChar(argin[0])) mexErrMsgTxt("argin error");
	if (mxGetString(argin[0],file,sizeof(file))!=0) mexErrMsgTxt("argin error");
	if (nargout>7) mexErrMsgTxt("argout error"); 
	
	nnav=ReadRinexNav(file,&data_b,&index,ionprm,dutc,comment,&type);
	if (nnav<0) {
		nsat=0;
		argout[0]=mxCreateCellArray(1,&nsat);
		argout[1]=mxCreateString("");
		argout[2]=mxCreateDoubleMatrix(0,0,mxREAL);
		argout[3]=mxCreateDoubleMatrix(0,0,mxREAL);
		argout[4]=mxCreateDoubleMatrix(0,0,mxREAL);
		argout[5]=mxCreateDoubleMatrix(0,0,mxREAL);
		argout[6]=mxCreateString("");
		return;
	}
	argout[0]=mxCreateCellArray(1,&nsat);
	for (i=1;i<=nsat;i++) {
		sprintf(sat,"%s%02d",type=='G'?"GLO":"GPS",i);
		mxSetCell(argout[0],i-1,mxCreateString(sat));
	}
	if ((p=strrchr(file,'\\'))!=NULL) p=p+1; else p=file;
	strncpy(rcv,p,4); rcv[4]='\0';
	argout[1]=mxCreateString(rcv);
	argout[2]=mxCreateDoubleMatrix(nnav,SIZE_NAV,mxREAL); data=mxGetPr(argout[2]);
	for (i=k=0;i<nnav;i++) for (j=0;j<SIZE_NAV;j++,k++) data[i+nnav*j]=data_b[k];
	argout[3]=mxCreateDoubleMatrix(nnav,1,mxREAL); COPY(index,nnav,1,mxGetPr(argout[3]));
	argout[4]=mxCreateDoubleMatrix(2,4,mxREAL); COPY(ionprm,2,4,mxGetPr(argout[4]));
	argout[5]=mxCreateDoubleMatrix(1,4,mxREAL); COPY(dutc,1,4,mxGetPr(argout[5]));
	argout[6]=mxCreateString(comment);
	free(data_b); free(index);
}
