/*------------------------------------------------------------------------------
% [system] : GpsTools
% [module] : read VMF1 grid file
% [func]   : read VMF1 grid file
% [argin]  : file = filepath
% [argout] : data  = data(i,j) lats(i) lons(j) grid data (nxm)
%            lats  = latitudes  (deg) (1xn)
%            lons  = longitudes (deg) (1xm)
% [note]   : reference :
%            J.Boehm et al., Troposphere mapping functions for GPS and very long
%            baseline interferometry from European Centre for Medium-Range
%            Weather Forecasts operational analysis data, J. Geoph. Res.,
%            Vol. 111, B02406, 2006
%            vmf1 grid file can be downloaded from :
%            http://www.hg.tuwien.ac.at/~ecmwf1/
% [version]: $Revision: 3 $ $Date: 06/07/08 1:01 $
%            Copyright(c) 2004-2008 by T.Takasu, all rights reserved
% [history]: 08/12/06  0.1  new
%-----------------------------------------------------------------------------*/
#include <stdio.h>
#include "mex.h"
#include "qtcmn.h"

/* read vmf grid file --------------------------------------------------------*/
extern int ReadVmfGrid(const char *file, int *siz, double **data, double **lats,
					   double **lons)
{
	FILE *fp;
	double val[12];
	int i=-1,j=0,k,n;
	char buff[256];
	
	if ((fp=fopen(file,"r"))==NULL) return 0;
	
	while (fgets(buff,sizeof(buff),fp)) {
		if (i<0) {
			if (sscanf(buff,"%lf %lf %lf %lf %lf %lf",
					   val,val+1,val+2,val+3,val+4,val+5)<6) continue;
			siz[0]=(int)((val[0]-val[1])/val[4]+0.5)+1;
			siz[1]=(int)((val[3]-val[2])/val[5]+0.5)+1;
			if (siz[0]<=0||siz[1]<=0) {
				fclose(fp);
				return 0;
			}
			*data=MAT(siz[0],siz[1]);
			*lats=MAT(siz[0],1);
			*lons=MAT(siz[1],1);
			if (!*data||!*lats||!*lons) {
				mexErrMsgTxt("memory allocation error");
				fclose(fp);
				return 0;
			}
			for (k=0;k<siz[0];k++) (*lats)[k]=val[0]-val[4]*k;
			for (k=0;k<siz[1];k++) (*lons)[k]=val[2]+val[5]*k;
			i=0;
			continue;
		}
		n=sscanf(buff,"%lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf",
				 val,val+1,val+2,val+3,val+4,val+5,val+6,val+7,val+8,val+9,
				 val+10,val+11);
		for (k=0;k<n&&j<siz[1];k++,j++) {
			(*data)[i+j*siz[0]]=val[k];
		}
		if (j>=siz[1]) {i++; j=0;}
		if (i>=siz[0]) break;
	}
	fclose(fp);
	return i>0;
}
/* mex interface -------------------------------------------------------------*/
extern void mexFunction(int nargout, mxArray *argout[], int nargin,
					    const mxArray *argin[])
{
	double *data,*lats,*lons;
	char file[1024];
	int i,j,siz[2]={0};
	
	if (nargin<1||!mxIsChar(argin[0])) mexErrMsgTxt("argin error");
	if (mxGetString(argin[0],file,sizeof(file))!=0) mexErrMsgTxt("argin error");
	if (nargout>3) mexErrMsgTxt("argout error"); 
	
	if (!ReadVmfGrid(file,siz,&data,&lats,&lons)) {
		if (nargout>0) argout[0]=mxCreateDoubleMatrix(0,0,mxREAL);
		if (nargout>1) argout[1]=mxCreateDoubleMatrix(0,0,mxREAL);
		if (nargout>2) argout[2]=mxCreateDoubleMatrix(0,0,mxREAL);
		return;
	}
	argout[0]=mxCreateDoubleMatrix(siz[0],siz[1],mxREAL);
	argout[1]=mxCreateDoubleMatrix(siz[0],1,mxREAL);
	argout[2]=mxCreateDoubleMatrix(siz[1],1,mxREAL);
	COPY(data,siz[0],siz[1],mxGetPr(argout[0]));
	COPY(lats,siz[0],1,mxGetPr(argout[1]));
	COPY(lons,siz[1],1,mxGetPr(argout[2]));
	
	free(data); free(lats); free(lons);
}
