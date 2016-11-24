/*------------------------------------------------------------------------------
% [system] : GpsTools
% [module] : mjd to calender date/time
% [func]   : convert mjd to calender date/time
% [argin]  : mjd  = modified julan date (day)
%           (sec) = seconds of the day
% [argout] : dt   = date/time [year,month,day,hour,min,sec]
% [note]   : valid after 1582/10/10
% [version]: $Revision: 4 $ $Date: 06/07/08 1:01 $
%            Copyright(c) 2004-2006 by T.Takasu, all rights reserved
% [history]: 03/02/23  0.1  new
%            04/12/17  0.2  fix bug of 14:19:60->14:20:00
%-----------------------------------------------------------------------------*/
#include "qtcmn.h"

/* mjd to calender date/time -------------------------------------------------*/
extern void MjdToCal(const double *mjd, const double *sec, double *dt)
{
	double a,q,b,c,e,f,g,h,Y,M,D;

	a=floor(mjd[0]+2468570.0);
	if (sec!=NULL) {a=a+floor(sec[0]/86400.0); q=MOD(sec[0],86400.0);}
	else q=(mjd[0]+2468570.0-a)*86400.0;
	b=floor(a/36524.25);
	c=a-floor(36524.25*b+0.75);
	e=floor((c+1.0)/365.25025);
	f=c-floor(365.25*e)+31.0;
	g=floor(f/30.59);
	h=floor(g/11.0);
	Y=100.0*(b-49.0)+e+h;
	M=g-12.0*h+2.0;
	D=f-floor(30.59*g);
	dt[0]=Y; dt[1]=M; dt[2]=D;
	dt[3]=floor(q/3600.0); q-=dt[3]*3600.0;
	dt[4]=floor(q/59.9999999999);
	dt[5]=q-dt[4]*59.9999999999;
}
