function prm=prm_gpssats(t)
%-------------------------------------------------------------------------------
% [system] : GpsTools
% [module] : read satellite parameters
% [func]   : read satellite parameters
% [argin]  : (t) = mjd-gpst (defalut: recent)
% [argout] : prm = satellite parameters
%                  prm(n,:)={name,type,mass,xsize,d0,y0,b0,z0}
% [note]   :
% [version]: $Revision: 12 $ $Date: 2008-11-25 10:02:15 +0900 (火, 25 11 2008) $
%            Copyright(c) 2004-2006 by T.Takasu, all rights reserved
% [history]: 04/10/28  0.1  new
%            05/06/13  0.2  add argin t
%            06/07/07  0.3  separate satellite parameter file
%-------------------------------------------------------------------------------
if nargin<1, t=99999999; end

[dirs,f]=fileparts(which(mfilename));
file=fullfile(dirs,'data','sats_params.txt');
if ~exist(file), prm={}; return, end
[f1,f2,f3,f4,f5,f6,f7,f8,f9]=...
    textread(file,'%s%f%f%f%f%f%f%f%f','delimiter',',','commentstyle','matlab');
prm=[f1,num2cell([f2,f3,f4,f5,f6,f7,f8])];
prm=prm(f9<=t,:);
[p,i]=unique(prm(:,1));
prm=prm(i,:);
