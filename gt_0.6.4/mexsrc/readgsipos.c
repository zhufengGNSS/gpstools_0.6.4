/*------------------------------------------------------------------------------
% [system] : GpsTools
% [module] : read GSI position estimation file
% [func]   : read GSI position estimation file
% [argin]  : file = filepath
% [argout] : epoch = epoch time [year,month,day,hour,min,sec]
%            time  = time vector relativ to epoch
%            poss  = position
%                    poss(n,1:3) : time(n) position [x,y,z] (m)
%            psigs = std. deviations
%                    psigs(n,1:3) : time(n) pos std.dev [dx,dy,dz] (m)
% [note]   :
% [version]: $Revision: 16 $ $Date: 2008-12-12 15:49:30 +0900 (金, 12 12 2008) $
%            Copyright(c) 2004-2008 by T.Takasu, all rights reserved
% [history]: 05/06/02  0.1  separated from readpos.m
%            08/11/30  0.2  suppress warning
%-----------------------------------------------------------------------------*/
#include <stdio.h>
#include <string.h>
#include "mex.h"
#include "qtcmn.h"

#define NMAX	512		/* max number of epoch in a file */

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
/* read rinex clock body -----------------------------------------------------*/
extern int ReadGsiPos(const char *file, double *epoch, double *time,
					  double *poss, double *psigs)
{
	FILE *fp;
	char buff[256];
	int i,start=0,np=0;
	double t[6],td,ts,td0=0.0,ts0=0.0;
	
	if ((fp=fopen(file,"rt"))==NULL) return -1;
	
	while (fgets(buff,sizeof(buff),fp)!=NULL) {
		if (!strncmp(buff,"+DATA",5)) {start=1; continue;}
		else if (!strncmp(buff,"-DATA",5)) break;
		if (!start) continue;
		if (sscanf(SubStr(buff,1,19),"%lf %lf %lf %lf:%lf:%lf",&t[0],&t[1],
		    &t[2],&t[3],&t[4],&t[5])!=6) continue;
		CalToMjd(t,&td,&ts,6);
		if (td0==0.0) {
			for (i=0;i<6;i++) epoch[i]=t[i];
			td0=td;
			ts0=ts;
		}
		time[np]=(td-td0)*86400.0+ts-ts0;
		for (i=0;i<3;i++) poss[np*3+i]=StrToNum(buff,21+i*18,17);
		np++;
	}
	fclose(fp);
	return np;
}
/* mex interface -------------------------------------------------------------*/
extern void mexFunction(int nargout, mxArray *argout[], int nargin,
					    const mxArray *argin[])
{
	static double time[NMAX],poss[NMAX*3],psigs[NMAX*3];
	double epoch[6];
	char file[256];
	int i,j,np;

	if (nargin<1||!mxIsChar(argin[0])) mexErrMsgTxt("argin error");
	if (mxGetString(argin[0],file,sizeof(file))!=0) mexErrMsgTxt("argin error");
	if (nargout>4) mexErrMsgTxt("argout error"); 
	
    np=ReadGsiPos(file,epoch,time,poss,psigs);
	if (np<=0) {
		if (nargout>0) argout[0]=mxCreateDoubleMatrix(0,0,mxREAL);
		if (nargout>1) argout[1]=mxCreateDoubleMatrix(0,0,mxREAL);
		if (nargout>2) argout[2]=mxCreateDoubleMatrix(0,0,mxREAL);
		if (nargout>3) argout[3]=mxCreateDoubleMatrix(0,0,mxREAL);
		return;
	}
	if (nargout>0) {
		argout[0]=mxCreateDoubleMatrix(1,6,mxREAL);
		COPY(epoch,1,6,mxGetPr(argout[0]));
	}
	if (nargout>1) {
		argout[1]=mxCreateDoubleMatrix(np,1,mxREAL);
		COPY(time,np,1,mxGetPr(argout[1]));
	}
	if (nargout>2) {
		argout[2]=mxCreateDoubleMatrix(np,3,mxREAL);
		for (i=0;i<3;i++)
		for (j=0;j<np;j++)
			*(mxGetPr(argout[2])+np*i+j)=poss[3*j+i];
	}
	if (nargout>3) {
		argout[3]=mxCreateDoubleMatrix(np,3,mxREAL);
		for (i=0;i<3;i++)
		for (j=0;j<np;j++)
			*(mxGetPr(argout[3])+np*i+j)=psigs[3*j+i];
	}
}
