function [latc,lonc]=geomtogeoc(latm,lonm)
%-------------------------------------------------------------------------------
% [system] : GpsTools
% [module] : geomagnetic to geocentric position
% [func]   : transform geomagnetic to geocentric position
% [argin]  : latm,lonm = geomagnetic latitude,longitude (rad)
% [argout] : latc,lonc = geocentric latitude,longitude (rad)
% [note]   : IGRF2000
% [version]: $Revision: 2 $ $Date: 06/07/08 14:11 $
%            Copyright(c) 2004-2006 by T.Takasu, all rights reserved
% [history]: 04/12/06  0.1  new
%-------------------------------------------------------------------------------
lat0=79.5; lon0=-71.6; % IGRF2000 magnetic north pole

rm=[cos(lonm)*cos(latm);sin(lonm)*cos(latm);sin(latm)];
rc=Rz(-lon0*pi/180)*Ry(lat0*pi/180-pi/2)*rm;
latc=asin(rc(3)); lonc=atan2(rc(2),rc(1));

% coordinate rotation matrix ---------------------------------------------------
function R=Ry(t), R=[cos(t),0,-sin(t);0,1,0;sin(t),0,cos(t)];
function R=Rz(t), R=[cos(t),sin(t),0;-sin(t),cos(t),0;0,0,1];
