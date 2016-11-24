function prm=prm_plotgpv_def
%-------------------------------------------------------------------------------
% [system] : GpsTools
% [module] : default parameters of plotgpv
% [func]   : default parameters of plotgpv
% [argin]  : none
% [argout] : prm = parameters struct
% [note]   :
% [version]: $Revision: 2 $ $Date: 06/07/14 20:58 $
% [history]: 04/03/24  0.1  new
%-------------------------------------------------------------------------------
prm.tstart  =[2004,1,1,0,0,0];
prm.tend    =[2004,1,1,0,0,0];
prm.tint    =6;
prm.ft      =0;
prm.src     ='rsm';
prm.type    ={'temp',''};
prm.height  =[nan,nan]; % nan : all
prm.ptype   ={'map','contl1','cbar'};
prm.view    =[0,90];
prm.vproj   =0;
prm.zlim    =[0.1,-inf,inf]; % inf:auto
prm.zpos    =[2,1,1];
prm.zfact   =0.1;
prm.axis    =1;
prm.box     =1;
prm.bcolor  ='w';
prm.cmap    =[];
prm.cbpos   =0;
prm.cbdata  =1;
prm.ccolor  ={'none','k'};
prm.clabel  ={0,0};
prm.cstep   ={[0,0,0],[0,0,0]};
prm.clevel  ={0,0};
prm.mcolor  =[0.5,0.5,0.5];
prm.mint    =2;
prm.fcolor  ='w';
prm.fparam  =[0.5,0.5,0.5];
prm.vcolor  ='b';
prm.vlen    =1;
prm.vint    =5;
prm.scolor  =[0.5,0.5,0.5];
prm.slen    =1000;
prm.sdir    =0;
prm.sint    =5;
prm.light   =[135,45];
prm.lstr    =[0.4,0.2];
prm.map.proj='lambert';
prm.map.cent=[135,40];
prm.map.base=[140,0];
prm.map.scale=1;
prm.map.gint=[10,0];
prm.map.color={'none','none','k','k'};
prm.dirs.gpv='';
