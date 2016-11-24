function prm=prm_gpsfiles_def
%-------------------------------------------------------------------------------
% [system] : GpsTools
% [module] : default parameters of gpsfiles
% [func]   : default parameters of gpsfiles
% [argin]  : none
% [argout] : prm = parameters struct
% [note]   :
% [version]: $Revision: 3 $ $Date: 06/07/18 14:10 $
% [history]: 04/07/02  0.1  new
%-------------------------------------------------------------------------------
prm.pos     =[0,0,600,400];
prm.possep  =0.2;
prm.dirs    ='';
prm.filter  ={'','',''};
prm.patt    ='*.*';
prm.sort    =1;
prm.ud      =0;
prm.compsufx=0;
% text viewer,image viewer,browser
prm.cmd1='';
prm.cmd2='';
prm.cmd3='c:\Program Files\Internet Explorer\iexplore.exe';
