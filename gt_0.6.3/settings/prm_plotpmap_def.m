function prm=prm_plotpmap_def
%-------------------------------------------------------------------------------
% [system] : GpsTools
% [module] : display/plot pwv estimations default parameters
% [func]   : display/plot pwv estimations default parameters
% [argin]  : none
% [argout] : prm = parameters
% [note]   :
% [version]: $Revision: 2 $ $Date: 06/07/14 20:58 $
% [history]: 05/11/24  0.1  new
%-------------------------------------------------------------------------------

prm.ep=[2004,1,1,0,0,0];
prm.ts=[0,1];
prm.tint=3600;
prm.rcvs={};
prm.area=[30,50,125,150];
prm.gint=0.1;
prm.intp=0;
prm.range=[0,0.5,80];
prm.showf=[0,0];
prm.cbpos=1;
prm.ccolor=[.9,.9,.9];
prm.pcolor='r';
prm.map.proj='miller';
prm.map.cent=[138,38];
prm.map.base=[138,0];
prm.map.scale=15;
prm.map.gint=[5,5];
prm.map.color={'none','w','k',[.5,.5,.5]};
prm.dirs.est='';
