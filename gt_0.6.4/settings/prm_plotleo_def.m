function prm=prm_plotleo_def
%-------------------------------------------------------------------------------
% [system] : GpsTools
% [module] : default parameters of plotleo
% [func]   : default parameters of plotleo
% [argin]  : none
% [argout] : prm = parameters struct
% [note]   :
% [version]: $Revision: 20 $ $Date: 2009-05-01 04:15:33 +0900 (金, 01 5 2009) $
% [history]: 04/03/24  0.1  new
%            08/12/15  0.2  add prm.refeph
%-------------------------------------------------------------------------------
prm.tstart  =[2004,1,1,0,0,0];
prm.tend    =[2004,1,1,0,0,0];
prm.tint    =300;
prm.tunit   =24;
prm.sats    ={};
prm.sat     ='';
prm.refeph  ='GRA00';
prm.range   =[-1,1];
prm.off     =[0,0,0];
prm.nsig    =2;
prm.fb      ='f';
prm.type    ={'poserr'};
prm.showf   =[0,0,0];
prm.ptype   =2;
prm.psize   =[1,5];
prm.pcolor  ='b';
prm.interp  =0;
prm.files   ={};
prm.dirs.est='';
prm.dirs.eph='';
prm.map.proj='miller';
prm.map.cent=[0,0];
prm.map.base=[0,0];
prm.map.scale=1;
prm.map.gint=[15,30];
prm.map.color={'none','none','k','k'};
