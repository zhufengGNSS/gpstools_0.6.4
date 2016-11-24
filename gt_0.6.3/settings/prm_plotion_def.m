function prm=prm_plotion_def
%-------------------------------------------------------------------------------
% [system] : GpsTools
% [module] : default parameters of plotion
% [func]   : default parameters of plotion
% [argin]  : none
% [argout] : prm = parameters struct
% [note]   :
% [version]: $Revision: 2 $ $Date: 06/07/14 20:58 $
% [history]: 04/03/24  0.1  new
%-------------------------------------------------------------------------------
prm.tstart  =[2004,1,1,0,0,0];
prm.tend    =[2004,1,1,0,0,0];
prm.tint    =3600;
prm.tunit   =24;
prm.fb      ='igs';
prm.ref     ='';
prm.ptype   =0;
prm.popt    =[0,0];
prm.ccolor  ='w';
prm.ecolor  ='w';
prm.fixlt   =0;
prm.range   =[0,2.5,100];
prm.light   =[135,45];
prm.lstr    =[0.5,0.2];
prm.map.proj='lambert';
prm.map.cent=[135,40];
prm.map.base=[140,0];
prm.map.scale=1;
prm.map.gint=[30,0];
prm.map.color={'none','none','k','k'};
prm.dirs.est='';
prm.dirs.ion='';
prm.dirs.nav='';
