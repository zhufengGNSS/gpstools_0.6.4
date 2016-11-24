/*------------------------------------------------------------------------------
% [system] : GpsTools
% [module] : relativity correction
% [func]   : relativity correction
% [argin]  : rsat = satellite postion(m) (eci)
%            vsat = satellite velocity(m/sec) (eci)
%            rrcv = station position(m) (eci)
%           (opts)= option flag (1=shapiro time delay correction) (default:1)
% [argout] : rels = relativity correction(m)
% [note]   :
% [version]: $Revision: 16 $ $Date: 2008-12-12 15:49:30 +0900 (金, 12 12 2008) $
%            Copyright(c) 2004-2006 by T.Takasu, all rights reserved
% [history]: 04/10/19   0.1  new
%            08/11/25   0.2  add argin opts (gt_0.6.4)
%-----------------------------------------------------------------------------*/
#include "qtcmn.h"

/* relativity correction -----------------------------------------------------*/
extern double RelCorr(const double *rsat, const double *vsat,
					  const double *rrcv, int opts)
{
	double rels,rs,rr,rrs,rsr[3];
	
	/* satellite clock correction */
	rels=2.0/M_C*DOT(rsat,vsat);
	
	if (opts==1) { /* shapiro time delay correction */
		rs=NORM(rsat); rr=NORM(rrcv);
		rsr[0]=rsat[0]-rrcv[0];
		rsr[1]=rsat[1]-rrcv[1];
		rsr[2]=rsat[2]-rrcv[2];
		rrs=NORM(rsr);
		rels+=2.0*GME/M_C/M_C*log((rs+rr+rrs)/(rs+rr-rrs));
	}
	return rels;
}
