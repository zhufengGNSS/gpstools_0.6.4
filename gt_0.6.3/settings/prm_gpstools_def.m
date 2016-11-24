function prm=prm_gpstools_def
%-------------------------------------------------------------------------------
% [system] : GpsTools
% [module] : default parameters of gpstools
% [func]   : default parameters of gpstools
% [argin]  : none
% [argout] : prm = parameters struct
% [note]   :
% [version]: $Revision: 3 $ $Date: 06/07/20 10:29 $
% [history]: 04/03/24  0.1  new
%-------------------------------------------------------------------------------
prm.pos     = [0,0,600,23];
prm.fonts   ={{'Arial',8},{'Arial',9},{'Courier New',9},{'Times New Roman',8}};
prm.maxtext =20000;
prm.editor  ='notepad.exe';
prm.browser ='c:\Program Files\Internet Explorer\iexplore.exe';

prm.genprod.tstart  =[2004,1,1,0,0,0];
prm.genprod.tend    =[2004,1,1,0,0,0];
prm.genprod.tint    =300;
prm.genprod.tunit   =24;
prm.genprod.type    ='eph';
prm.genprod.fb      ='fb';
prm.genprod.opt     ='';
prm.genprod.list    ={{},{}};
prm.genprod.idir    ='';
prm.genprod.odir    ='';
prm.genprod.prmfile ='';

prm.geneph.tstart  =[2004,1,1,0,0,0];
prm.geneph.tend    =[2004,1,1,0,0,0];
prm.geneph.tint    =30;
prm.geneph.tunit   =24;
prm.geneph.eph     ='igs';
prm.geneph.erp     ='igs';
prm.geneph.ephdir  ='';
prm.geneph.erpdir  ='';
prm.geneph.odir    ='';

prm.genpwv.tstart  =[2004,1,1,0,0,0];
prm.genpwv.tend    =[2004,1,1,0,0,0];
prm.genpwv.tint    =30;
prm.genpwv.tunit   =24;
prm.genpwv.rcvs    ={};
prm.genpwv.estdir  ='';
prm.genpwv.metdir  ='';
prm.genpwv.odir    ='';

prm.genoload.time  =[2004,1,1,0,0,0];
prm.genoload.rcvs  ={};
prm.genoload.psrc  ='approx';
prm.genoload.tunit =24;
prm.genoload.gdir  ='%P\..\GOTIC2';
prm.genoload.pdir  ='';
prm.genoload.file  ='output.blq';

prm.multp.tstart   =[2004,1,1,0,0,0];
prm.multp.tend     =[2004,1,1,0,0,0];
prm.multp.tint     =30;
prm.multp.tunit    =24;
prm.multp.fb       ='b';
prm.multp.rcvs     ={};
prm.multp.nmax     =12;
prm.multp.estdir   ='';
prm.multp.navdir   ='';
prm.multp.odir     ='';
