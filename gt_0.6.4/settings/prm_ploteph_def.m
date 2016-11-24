function prm=prm_ploteph_def
%-------------------------------------------------------------------------------
% [system] : GpsTools
% [module] : default parameters of ploteph
% [func]   : default parameters of ploteph
% [argin]  : none
% [argout] : prm = parameters struct
% [note]   :
% [version]: $Revision: 12 $ $Date: 2008-11-25 10:02:15 +0900 (火, 25 11 2008) $
% [history]: 04/03/24  0.1  new
%-------------------------------------------------------------------------------
prm.tstart  =[2004,1,1,0,0,0];
prm.tend    =[2004,1,1,0,0,0];
prm.tint    =300;
prm.tunit   =24;
prm.sats    ={};
prm.sat     ='';
prm.range   =[-0.5,0.5];
prm.nsig    =0;
prm.fb      ='ephf';
prm.type    ={'pos3d','posrac'};
prm.ref     ='igs';
prm.interp  =0;
prm.showf   =[1,1];
prm.ptype   =2;
prm.psize   =[1,5];
prm.files   ={};
prm.dirs.est='';
prm.dirs.eph='';
prm.dirs.clk='';
prm.dirs.erp='';
