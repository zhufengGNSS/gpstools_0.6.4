/*------------------------------------------------------------------------------
% [system] : GpsTools
% [module] : relativity correction
% [func]   : relativity correction
% [argin]  : rsat = satellite postion(m) (eci)
%            vsat = satellite velocity(m/sec) (eci)
%            rrcv = station position(m) (eci)
% [argout] : rels = relativity correction(m)
% [note]   :
% [version]: $Revision: 4 $ $Date: 06/07/08 1:01 $
%            Copyright(c) 2004-2006 by T.Takasu, all rights reserved
% [history]: 04/10/19   0.1  new
%-----------------------------------------------------------------------------*/
#include "qtcmn.h"

/* relativity correction -----------------------------------------------------*/
extern double RelCorr(const double *rsat, const double *vsat,
					  const double *rrcv)
{
	double rels,rs,rr,rrs,rsr[3];
	
	/* satellite clock correction */
	rels=2.0/M_C*DOT(rsat,vsat);
	
	/* signal propagation correction */
	rs=NORM(rsat); rr=NORM(rrcv);
	rsr[0]=rsat[0]-rrcv[0];
	rsr[1]=rsat[1]-rrcv[1];
	rsr[2]=rsat[2]-rrcv[2];
	rrs=NORM(rsr);
	rels+=2.0*GME/M_C/M_C*log((rs+rr+rrs)/(rs+rr-rrs));
	return rels;
}
