function prm=prm_plottrop_def
%-------------------------------------------------------------------------------
% [system] : GpsTools
% [module] : default parameters of plottrop
% [func]   : default parameters of plottrop
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
prm.range   =[-0.05,0.05];
prm.nsig    =2;
prm.fb      ='zpdf';
prm.type    ={'ztd'};
prm.ref     ='igs';
prm.map.proj='miller';
prm.map.cent=[0,0];
prm.map.base=[0,0];
prm.map.scale=1;
prm.map.gint=[15,30];
prm.map.color={'none','none','k','k'};
prm.showf   =[0,1];
prm.ptype   =2;
prm.psize   =[1,5];
prm.files   ={};
prm.dirs.est ='';
prm.dirs.trop='';
