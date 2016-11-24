function prm=prm_plotstats_def
%-------------------------------------------------------------------------------
% [system] : GpsTools
% [module] : default parameters of plotstats
% [func]   : default parameters of plotstats
% [argin]  : none
% [argout] : prm = parameters struct
% [note]   :
% [version]: $Revision: 3 $ $Date: 06/07/14 20:58 $
% [history]: 04/03/24  0.1  new
%-------------------------------------------------------------------------------
prm.tstart  =[2004,1,1,0,0,0];
prm.tend    =[2004,1,1,0,0,0];
prm.tint    =300;
prm.tunit   =24;
prm.srcnav  ='rinex';
prm.fb      ='resf';
prm.type    ='res';
prm.sats    ={};
prm.sat     ='ALL';
prm.rcvs    ={};
prm.rcv     ='ALL';
prm.execl   =0;
prm.ptype   =1;
prm.psize   =5;
prm.mesh    =[5,5];
prm.range   =[-0.05,0.05];
prm.showf   =[1,1,1,1];
prm.vangle  =[45,30];
prm.files   ={};
prm.dirs.est='';
prm.dirs.nav='';
