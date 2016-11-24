/*------------------------------------------------------------------------------
% [system] : GpsTools
% [module] : read rinex observation data
% [func]   : read rinex observation data file
% [argin]  : file    = file path
% [argout] : epoch   = first obs. epoch [year,month,day,hour,min,sec]
%            time    = time vector relative to epoch (sec)
%            types   = observation data types
%            units   = observation data units
%            sats    = satellites list
%            rcv     = receiving station
%            data    = observation data
%                      data(n,m) : time(n),index(n),types(m) data
%            index   = satellite index
%            rcvpos  = approx station position (m) [x;y;z]
%            antdel  = antenna delta (m) [up;east;north]
%            anttype = antenna model name
%            rcvtype = receiver model name
%            comment = rinex header lines separated by ';'
% [note]   :
% [version]: $Revision: 12 $ $Date: 06/07/08 1:01 $
%            Copyright(c) 2004-2006 by T.Takasu, all rights reserved
% [history]: 04/11/16  0.1  new
%            05/08/15  0.2  recognize 0 as missing observation
%            05/09/30  0.3  data valid even if epoch flag != 0
%            05/10/12  0.4  fix bug to crash if header comment size overflows
%-----------------------------------------------------------------------------*/
#include <stdio.h>
#include <string.h>
#include "mex.h"
#include "qtcmn.h"

#define	COMSIZMAX		4096		/* comment max size */

extern void CalToMjd(const double *dt, double *mjd, double *sec, int len_dt);

/* string to number ----------------------------------------------------------*/
static double StrToNum(const char *str, int pos, int len)
{
	double value;
	char buff[64];
	if (strlen(str)<pos+len) return mxGetNaN();
	strncpy(buff,str+pos,len); buff[len]='\0';
	if (sscanf(buff,"%lf",&value)==1) return value;
	return mxGetNaN();
}
/* read rinex header ---------------------------------------------------------*/
static int ReadHead(FILE *fp, char *type, double *epoch, double *rcvpos,
				    double *antdel,char *anttype, char *rcvtype, char **types,
				    char **units, char *comment)
{
	static char *rtypes[]=
		{"L1","L2","C1","P1","P2","D1","D2","T1","T2","S1","S2"};
	static char *runits[]=
		{"cycle","cycle","m","m","m","Hz","Hz","cycle","cycle","",""};
	int i,j,nobs=0;
	char buff[256],*label,*p=comment;

	while (fgets(buff,sizeof(buff),fp)!=NULL) {
		buff[255]='\0';
		if (p+strlen(buff)+1<comment+COMSIZMAX) {
			strcpy(p,buff); p+=strlen(p)-1; strcpy(p++,";");
		}
		label=buff+60;
		
		if (strstr(label,"RINEX VERSION / TYPE")!=NULL) {
			*type=buff[20];
		}
		else if (strstr(label,"# / TYPES OF OBSERV")!=NULL) {
			if (sscanf(buff,"%d",&nobs)!=1) continue;
			for (i=0;i<nobs;i++) {
				if (i!=0&&i%9==0) fgets(buff,sizeof(buff),fp);
				for (j=0,types[i]=units[i]="";j<11;j++) {
					if (strncmp(buff+10+(i%9)*6,rtypes[j],2)!=0) continue;
					types[i]=rtypes[j];
					units[i]=runits[j];
				}
			}
		}
		else if (strstr(label,"APPROX POSITION XYZ")!=NULL) {
			for (i=0;i<3;i++) rcvpos[i]=StrToNum(buff,i*14,14);
		}
		else if (strstr(label,"ANTENNA: DELTA H/E/N")!=NULL) {
			for (i=0;i<3;i++) antdel[i]=StrToNum(buff,i*14,14);
		}
		else if (strstr(label,"ANT # / TYPE")!=NULL) {
			strncpy(anttype,buff+20,20);
			for (i=19;i>=0;i--) if (anttype[i]!=' ') break;
			anttype[i+1]='\0';
		}
		else if (strstr(label,"REC # / TYPE")!=NULL) {
			strncpy(rcvtype,buff+20,20);
			for (i=19;i>=0;i--) if (rcvtype[i]!=' ') break;
			rcvtype[i+1]='\0';
		}
		else if (strstr(label,"TIME OF FIRST OBS")!=NULL) {
			if (sscanf(buff,"%lf %lf %lf %lf %lf %lf",&epoch[0],&epoch[1],
					   &epoch[2],&epoch[3],&epoch[4],&epoch[5])!=6) return 0;
		}
		else if (strstr(label,"END OF HEADER")!=NULL) break;
	}
	return nobs;
}
/* read rinex body -----------------------------------------------------------*/
static int ReadData(FILE *fp, double td, double ts, int nobs, double **time,
				    double **data, double **index)
{
	double epoch[6],t,tdd,tss,value;
	int i,j,flag,nsat,nz=0,size=4096,sat[64];
	char buff[256],*id;
	
	*time =(double *)malloc(size*sizeof(double));
	*index=(double *)malloc(size*sizeof(double));
	*data =(double *)malloc(size*nobs*sizeof(double));
	if (*time==NULL||*index==NULL||*data==NULL) return 0;
	
	while (fgets(buff,sizeof(buff),fp)!=NULL) {
		if (sscanf(buff,"%lf %lf %lf %lf %lf %lf",&epoch[0],&epoch[1],
				   &epoch[2],&epoch[3],&epoch[4],&epoch[5])!=6) continue;
		epoch[0]+=epoch[0]<80.0?2000.0:1900.0;
		CalToMjd(epoch,&tdd,&tss,6);
		t=(tdd-td)*86400.0+tss-ts;
		
		flag=StrToNum(buff,26,3);
		nsat=StrToNum(buff,29,3);
		for (i=0;i<nsat;i++) {
			if (i!=0&&i%12==0) fgets(buff,sizeof(buff),fp);
			id=buff+32+(i%12)*3;
			if (*id!=' '&&*id!='G') sat[i]=0;
			else sat[i]=StrToNum(buff,33+(i%12)*3,2);
		}
		for (i=0;i<nsat;i++) {
			if (nz>=size) {
				size*=2;
				*time =(double *)realloc(*time,size*sizeof(double));
				*index=(double *)realloc(*index,size*sizeof(double));
				*data =(double *)realloc(*data,size*nobs*sizeof(double));
				if (*time==NULL||*index==NULL||*data==NULL) return 0;
			}
			(*time)[nz]=t;
			(*index)[nz]=(double)sat[i];
			for (j=0;j<nobs;j++) {
				if (j%5==0) fgets(buff,sizeof(buff),fp);
				value=StrToNum(buff,(j%5)*16,14);
				if (ABS(value)<0.001) value=mxGetNaN();
				(*data)[nz*nobs+j]=value;
			}
			if (sat[i]!=0) nz++;
		}
	}
	return nz;
}
/* read rinex obs data -------------------------------------------------------*/
extern int ReadRinexObs(const char *file, double *epoch, double **time,
					    char **types, char **units, double **data,
					    double **index, double *rcvpos, double *antdel,
					    char *anttype, char *rcvtype, char *comment,int *nobs)
{
	FILE *fp;
	char type;
	double td,ts;
	int nz;

	if ((fp=fopen(file,"rt"))==NULL) return -1;
	*nobs=ReadHead(fp,&type,epoch,rcvpos,antdel,anttype,rcvtype,types,units,comment);
	if (*nobs<=0||type!='O') {
		fclose(fp);
		return -1;
	}
	CalToMjd(epoch,&td,&ts,6);
	nz=ReadData(fp,td,ts,*nobs,time,data,index);
	fclose(fp);
	return nz;
}
/* mex interface -------------------------------------------------------------*/
extern void mexFunction(int nargout, mxArray *argout[], int nargin,
					    const mxArray *argin[])
{
	double epoch[6],*time,rcvpos[3],antdel[3],*data,*index,*time_b,*data_b,*index_b;
	char *p,file[256],*types[32],*units[32],comment[COMSIZMAX]="",sat[16],rcv[16];
	char anttype[32]="",rcvtype[32]="";
	int i,j,k,nz,nobs,nsat=31;
	
	if (nargin<1||!mxIsChar(argin[0])) mexErrMsgTxt("argin error");
	if (mxGetString(argin[0],file,sizeof(file))!=0) mexErrMsgTxt("argin error");
	if (nargout<8||13<nargout) mexErrMsgTxt("argout error"); 
	
	nz=ReadRinexObs(file,epoch,&time_b,types,units,&data_b,&index_b,rcvpos,
					antdel,anttype,rcvtype,comment,&nobs);
	if (nz<0) {
		argout[0]=mxCreateDoubleMatrix(0,0,mxREAL); 
		argout[1]=mxCreateDoubleMatrix(0,0,mxREAL);
		nobs=0;
		argout[2]=mxCreateCellArray(1,&nobs);
		argout[3]=mxCreateCellArray(1,&nobs);
		argout[4]=mxCreateCellArray(1,&nobs);
		argout[5]=mxCreateString("");
		argout[6]=mxCreateDoubleMatrix(0,0,mxREAL);
		argout[7]=mxCreateDoubleMatrix(0,0,mxREAL);
		if (nargout>8) argout[8]=mxCreateDoubleMatrix(0,0,mxREAL);
		if (nargout>9) argout[9]=mxCreateDoubleMatrix(0,0,mxREAL);
		if (nargout>10) argout[10]=mxCreateString("");
		if (nargout>11) argout[11]=mxCreateString("");
		if (nargout>12) argout[12]=mxCreateString("");
		return;
	}
	argout[0]=mxCreateDoubleMatrix(1,6,mxREAL); COPY(epoch,1,6,mxGetPr(argout[0]));
	argout[1]=mxCreateDoubleMatrix(nz,1,mxREAL); time=mxGetPr(argout[1]);
	COPY(time_b,nz,1,time);
	argout[2]=mxCreateCellArray(1,&nobs);
	argout[3]=mxCreateCellArray(1,&nobs);
	for (i=0;i<nobs;i++) {
		mxSetCell(argout[2],i,mxCreateString(types[i]));
		mxSetCell(argout[3],i,mxCreateString(units[i]));
	}
	argout[4]=mxCreateCellArray(1,&nsat);
	for (i=1;i<32;i++) {
		sprintf(sat,"GPS%02d",i);
		mxSetCell(argout[4],i-1,mxCreateString(sat));
	}
	if ((p=strrchr(file,'\\'))!=NULL) p=p+1; else p=file;
	strncpy(rcv,p,4); rcv[4]='\0';
	argout[5]=mxCreateString(rcv);
	argout[6]=mxCreateDoubleMatrix(nz,nobs,mxREAL); data=mxGetPr(argout[6]);
	for (i=k=0;i<nz;i++) for (j=0;j<nobs;j++,k++) data[i+nz*j]=data_b[k];
	argout[7]=mxCreateDoubleMatrix(nz,1,mxREAL); index=mxGetPr(argout[7]);
	COPY(index_b,nz,1,index);
	if (nargout>8) {
		argout[8]=mxCreateDoubleMatrix(3,1,mxREAL);
		COPY(rcvpos,3,1,mxGetPr(argout[8]));
	}
	if (nargout>9) {
		argout[9]=mxCreateDoubleMatrix(3,1,mxREAL);
		COPY(antdel,3,1,mxGetPr(argout[9]));
	}
	if (nargout>10) argout[10]=mxCreateString(anttype);
	if (nargout>11) argout[11]=mxCreateString(rcvtype);
	if (nargout>12) argout[12]=mxCreateString(comment);
	
	free(time_b); free(data_b); free(index_b);
}
