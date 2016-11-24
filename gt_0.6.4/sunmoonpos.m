function [rsun,rmoon]=sunmoonpos(tutc)
%-------------------------------------------------------------------------------
% [system] : GpsTools
% [module] : sun/moon position
% [func]   : get solar/lunar positions
% [argin]  : tutc  = time (mjd-utc)
% [argout] : rsun  = solar position (m) (eci)
%            rmoon = lunar position (m) (eci)
% [note]   :
% [version]: $Revision: 12 $ $Date: 2008-11-25 10:02:15 +0900 (火, 25 11 2008) $
%            Copyright(c) 2004-2006 by T.Takasu, all rights reserved
% [history]: 04/11/19  0.1  new
%-------------------------------------------------------------------------------
global p_ ephpdir
if isempty(ephpdir)
    [dirs,f]=fileparts(which(mfilename));
    ephpdir=fullfile(dirs,'data');
end
r=ephpl(tutc,[3,10,11])*1E3; % read eath-moon mass center,moon and sun pos.
rmoon=r(:,2);
rsun=r(:,3)+rmoon/82.30056-r(:,1);
