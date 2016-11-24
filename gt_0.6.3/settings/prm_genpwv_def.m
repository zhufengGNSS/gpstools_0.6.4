function prm=prm_genpwv_def
%-------------------------------------------------------------------------------
% [system] : GpsTools
% [module] : default parameters of generate pwv
% [func]   : default parameters of generate pwv
% [argin]  : none
% [argout] : prm = parameters struct
% [note]   :
% [version]: $Revision: 2 $ $Date: 06/07/14 20:58 $
% [history]: 04/07/25  0.1  new
%-------------------------------------------------------------------------------
prm.time    =[2004,1,1,0,0,0];
prm.tint    =3600;
prm.span    =1;
prm.rcvs    ={};
prm.tunit   =24;
prm.dirs.trop='';
prm.dirs.met ='';
prm.outdir   ='';
