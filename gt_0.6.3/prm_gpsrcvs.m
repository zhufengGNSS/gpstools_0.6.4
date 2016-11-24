function prm=prm_gpsrcvs
%------------------------------------------------------------------------------
% [system] : GpsTools
% [module] : read receiver definitions
% [func]   : read receiver definitions
% [argin]  : none
% [argout] : prm : receiver definitions
%                  prm(n,:)={name,category,fullname,comment}
% [note]   :
% [version]: $Revision: 9 $ $Date: 06/07/08 14:17 $
%            Copyright(c) 2004-2006 by T.Takasu, all rights reserved
% [history]: 04/11/03  0.1  new
%            06/07/07  0.2  separate receiver parameter file
%------------------------------------------------------------------------------
[dirs,f]=fileparts(which(mfilename));
file=fullfile(dirs,'data','rcvs_params.txt');
if ~exist(file), prm={}; return, end
[f1,f2,f3,f4]=textread(file,'%s%s%s%s%*[^\n]','delimiter',',',...
                       'commentstyle','matlab');
prm=[f1,f2,f3,f4];
