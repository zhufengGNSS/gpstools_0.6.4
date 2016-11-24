function prm=prm_plotclk_def
%-------------------------------------------------------------------------------
% [system] : GpsTools
% [module] : default parameters of plotclk
% [func]   : default parameters of plotclk
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
prm.sr      ='sat';
prm.sats    ={};
prm.sat     ='';
prm.rcvs    ={};
prm.rcv     ='';
prm.range   =[-1,1];
prm.nsig    =2;
prm.fb      ='clkf';
prm.type    ={'clkbias','clkerr'};
prm.ref     ='igs';
prm.refclk  ='';
prm.showf   =[0,1,0];
prm.adjc    =0;
prm.interp  =0;
prm.ptype   =2;
prm.psize   =[1,5];
prm.sec     =1;
prm.dirs.est='';
prm.dirs.clk='';
