/*------------------------------------------------------------------------------
% [system] : GpsTools
% [module] : calender date/time to mjd
% [func]   : convert calender date/time to mjd
% [argin]  : dt  = date/time [year,month,day(,hour,min,sec)]
% [argout] : mjd = modefied julan date (day)
%            (sec) = sec
% [note]   : valid after 1582/10/10
% [version]: $Revision: 12 $ $Date: 2008-11-25 10:02:15 +0900 (火, 25 11 2008) $
%            Copyright(c) 2004-2006 by T.Takasu, all rights reserved
% [history]: 04/02/23  0.1  new
%-----------------------------------------------------------------------------*/
#include "qtcmn.h"

/* calender date/time to mjd -------------------------------------------------*/
extern void CalToMjd(const double *dt, double *mjd, double *sec, int len_dt)
{
	double dt0=dt[0], dt1=dt[1];

	if (dt1<=2.0) {dt0=dt[0]-1.0; dt1=dt1+12.0;}
	mjd[0]=365.0*dt0-679004.0+floor(dt0/400.0)-floor(dt0/100.0)+
		   floor(dt0/4.0)+floor(30.6001*(dt1+1.0))+dt[2];
	if (len_dt>=6) {
	    if (sec!=NULL) sec[0]=dt[3]*3600.0+dt[4]*60.0+dt[5];
	    else mjd[0]=mjd[0]+dt[3]/24.0+dt[4]/1440.0+dt[5]/86400.0;
	}
}
