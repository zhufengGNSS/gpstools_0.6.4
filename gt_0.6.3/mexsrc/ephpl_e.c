/*------------------------------------------------------------------------------
% [system] : GpsTools
% [module] : solar/planetary ephemeris
% [func]   : read solar/planetary postions from de405 ephemeris
% [argin]  : t = date/time(mjd-utc)
%            pl= solar/planetary numbers
%                (1:mercury,2:venus,3:earth-moon barycenter,4:mars,5:jupiter,
%                 6:saturn,7:uranus,8:neptune,9:pluto,10:moon,11:sun)
%            <global>
%            utc_tai = utc-tai(sec)
%            ephpdir = ephemris data directory
% [argout] : r = solar/planetary positions(km)
% [note]   : coordinate:barycenter(all except for moon),eci(icrf)(moon)
% [version]: $Revision: 4 $ $Date: 06/07/08 1:01 $
%            Copyright(c) 2004-2006 by T.Takasu, all rights reserved
% [history]: 04/03/13  0.1  new
%-----------------------------------------------------------------------------*/
#include <stdio.h>
#include "qtcmn.h"

#define EPHFILE		"ephem2000.405"		/* ephemeris file name */
#define RSIZE		1018				/* size of record */
#define RLEN		(RSIZE*8)			/* length of record(byte) */

typedef struct {						/* de header type */
	long	coef[12][3];
	long	deno;
	long	lib[3];
} DeHeadType;

/* global variables ----------------------------------------------------------*/
extern double *utc_tai;
char ephpdir[256]="";

/* local variables -----------------------------------------------------------*/
static char	Labels[3*84+400*6];
static double Times[3];
static long Nconst;
static double Const[2];
static DeHeadType De;					/* de header cache */
static double R[RSIZE];					/* de record cache */

/* read de ephemeris record --------------------------------------------------*/
static int ReadRecord(double t)
{
	FILE *fp;
	long pos;
	char file[256];

	/* make ephemeris file path */
	if (ephpdir[0]=='\0') strcpy(file,EPHFILE);
	else sprintf(file,"%s%s%s",ephpdir,FILESEP,EPHFILE);
	
	if ((fp=fopen(file,"rb"))==NULL) return 0;

	/* read some informations */
	fread(Labels, sizeof(Labels),1,fp);
	fread(Times,  sizeof(Times), 1,fp);
	fread(&Nconst,sizeof(Nconst),1,fp);
	fread(Const,  sizeof(Const), 1,fp);
	
	/* read de ephemeris header */
	if (fread(&De,sizeof(De),1,fp)!=1||fseek(fp,RLEN*2,SEEK_SET)!=0) {
		fclose(fp);
		return 0;
	}
	/* read de ephemeris record */
	if (fread(R,sizeof(double),RSIZE,fp)!=RSIZE) {
		fclose(fp);
		return 0;
	}
	pos=(long)((t-R[0])/(R[1]-R[0]))+2;
	if (pos<2||fseek(fp,RLEN*pos,SEEK_SET)!=0) {
		fclose(fp);
		return 0;
	}
	fread(R,sizeof(double),RSIZE,fp);
	fclose(fp);
	return R[0]<=t&&t<R[1];
}
/* solar/planetary ephemeris -------------------------------------------------*/
extern int EphPl(const double *t, const double *pl, int np, double *r)
{
	long i,j,k,n,nc,ng,pos;
	double tt,tspan,t0,tcc,tc[1024];
	
	/* terestrial time in jd */
	tt=t[0]+(32.184-utc_tai[0])/86400.0+2400000.5;

	/* read de ephemeris record */
	if (R[0]<=0.0||tt<R[0]||R[1]<=tt) if (!ReadRecord(tt)) return 0;
	
	for (i=0;i<np;i++) {
		if ((int)pl[i]<1||11<(int)pl[i]) continue;
		nc=De.coef[(int)pl[i]-1][1];
		ng=De.coef[(int)pl[i]-1][2];
		
		tspan=(R[1]-R[0])/(double)ng;
		n=(long)(tt-R[0])/tspan;
		t0=R[0]+n*tspan;
		pos=De.coef[(int)pl[i]-1][0]+nc*3*n-1;
		
		/* chebyshev polynomial */
		tcc=2.0*(tt-t0)/tspan-1.0; tc[0]=1.0; tc[1]=tcc;
		for (j=2;j<nc;j++) tc[j]=2.0*tcc*tc[j-1]-tc[j-2];
		for (j=0;j<3;j++) {
			r[i*3+j]=0.0;
			for (k=0;k<nc;k++) r[i*3+j]+=R[pos+nc*j+k]*tc[k];
		}
	}
	return 1;
}
