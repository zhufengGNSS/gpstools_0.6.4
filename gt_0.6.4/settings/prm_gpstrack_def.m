function prm=prm_gpstrack_def
%-------------------------------------------------------------------------------
% [system] : GpsTools
% [module] : default parameters of gpstrack
% [func]   : default parameters of gpstrack
% [argin]  : none
% [argout] : prm = parameters struct
% [note]   :
% [version]: $Revision: 12 $ $Date: 2008-11-25 10:02:15 +0900 (火, 25 11 2008) $
% [history]: 04/03/24  0.1  new
%-------------------------------------------------------------------------------
prm.tstart  =[2004,1,1,0,0,0];
prm.tend    =[2004,1,1,0,0,0];
prm.tint    =900;
prm.sats    ={};
prm.sat     ='ALL';
prm.rcvs    ={};
prm.rcv     ='ALL';
prm.type    ={};
prm.showf   =[0,1,0];
prm.ptype   =2;
prm.elmin   =10;
prm.dirs.nav='';
prm.dirs.pos='';
prm.src.nav ='brdc';
prm.rcolor  ='r';
prm.rfont   ='';
prm.psize   =[1,5];
prm.map.proj='miller';
prm.map.cent=[0,0];
prm.map.base=[0,0];
prm.map.scale=1;
prm.map.gint=[15,30];
prm.map.color={'w','none',[0.5,0.5,0.5],[0.5,0.5,0.5]};
