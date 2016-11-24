function prm=prm_plotpos_def
%-------------------------------------------------------------------------------
% [system] : GpsTools
% [module] : default parameters of plotpos
% [func]   : default parameters of plotpos
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
prm.rcvs    ={};
prm.rcv     ='';
prm.range   =[-0.2,0.2];
prm.nsig    =2;
prm.fb      ='posfb';
prm.type    ='post';
prm.ref     ='estm';
prm.detr    =0;
prm.rfont   ='';
prm.fitn    =2;
prm.map.proj='miller';
prm.map.cent=[0,0];
prm.map.base=[0,0];
prm.map.scale=1;
prm.map.gint=[15,30];
prm.map.color={'w','none',[0.5,0.5,0.5],[0.5,0.5,0.5]};
prm.est     =0; % 0:all estimation, 1:final estimation only
prm.showf   =[1,1,0,0];
prm.adj     ='';
prm.sigmax  =0;
prm.srday   =[0,0,0];
prm.srprm   =[10,10];
prm.bpfilt  =[0,0];
prm.ptype   =1;
prm.psize   =[1,5];
prm.trprm   ='';
prm.files   ={};
prm.dirs.est='';
prm.dirs.pos='';
