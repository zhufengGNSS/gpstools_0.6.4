function h=geodh(gpos,model)
%-------------------------------------------------------------------------------
% [system] : GpsTools
% [module] : geodetic height
% [func]   : geodetic height
% [argin]  : gpos = geodetic position [lat,lon,h]
%                gpos(1) = latitude (deg)
%                gpos(2) = longitude (deg)
%                gpos(3) = ellipsoid height (m)
%           (model) = geoid model (default:'egm96')
%                'egm96' : egm96
%                'gsi'   : gsi2000
% [argout] : h    = geoidetic height (m)
% [note]   :
% [version]: $Revision: 5 $ $Date: 06/07/08 14:11 $
%            Copyright(c) 2004-2006 by T.Takasu, all rights reserved
% [history]: 05/06/15   0.1  new
%            05/07/09   0.2  substract geoid height from elipsoid height
%-------------------------------------------------------------------------------
persistent gprm data file
if nargin<2, model='egm96'; end
switch model, case 'egm96', f='geoid_egm96'; case 'gsi', f='geoid_gsi2000'; end
if isempty(data)|~strcmp(f,file)
    [dirs,e]=fileparts(which(mfilename));
    fs=fullfile(dirs,'data',f);
    load(fs); file=f;
end
lats=gprm.lat1:-gprm.dy:gprm.lat2;
lons=gprm.lon1:gprm.dx:gprm.lon2;
h=gpos(3)-interp2(lons,lats,double(data),gpos(2),gpos(1));
