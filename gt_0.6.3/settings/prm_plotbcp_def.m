function prm=prm_plotbcp_def
%-------------------------------------------------------------------------------
% [system] : GpsTools
% [module] : default parameters of plotbcp
% [func]   : default parameters of plotbcp
% [argin]  : none
% [argout] : prm = parameters struct
% [note]   :
% [version]: $Revision: 2 $ $Date: 06/07/14 20:58 $
% [history]: 04/03/24  0.1  new
%-------------------------------------------------------------------------------
prm.tstart  =[2004,1,1,0,0,0];
prm.tend    =[2004,1,1,0,0,0];
prm.tint    =300;
prm.tunit   =24;
prm.rcvs    ={};
prm.rcv     ='';
prm.range   =[-0.05,0.05];
prm.nsig    =0;
prm.type    ='';
prm.fb      ='bcpf';
prm.showf   =[1,1,0,0];
prm.ptype   =1;
prm.psize   =[1,5];
prm.dirs.est='';
