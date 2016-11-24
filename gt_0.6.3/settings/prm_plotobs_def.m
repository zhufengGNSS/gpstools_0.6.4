function prm=prm_plotobs_def
%-------------------------------------------------------------------------------
% [system] : GpsTools
% [module] : default parameters of plotobs
% [func]   : default parameters of plotobs
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
prm.type    ={'Data Availability','L1','L2'};
prm.sats    ={};
prm.sat     ='';
prm.rcvs    ={};
prm.rcv     ='';
prm.src     ='rinex';
prm.navsrc  ='brdc';
prm.nobias  =0;
prm.plotf   =3;
prm.showf   =[0,0,0];
prm.ptype   =1;
prm.color   ={[0.5,0.5,0.5],'b','r'};
prm.mgap    =60;
prm.range   =[-inf,inf];
prm.teqcopt ='+qc +set_mask 0 +ssv +sym +l -plot';
prm.dirs.obs='';
prm.dirs.obc='';
prm.dirs.nav='';
