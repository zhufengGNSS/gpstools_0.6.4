function trop=trop_mhop(t,azel,gpos,met_prm)
%-------------------------------------------------------------------------------
% [system] : GpsTools
% [module] : tropospheric model - modefied Hopfield model
% [func]   : calculate tropospheric delay by modefied Hopfield model
% [argin]  : t = time (mjd-utc) (no use)
%            azel = satellite azimath/elevation angle (rad) [az,el]
%            gpos = latitude/longitude/height (deg,m) [lat,lon,h]
%           (met_prm) = meteological parameters [pres,temp,humi]
%                 pres = pressure (hPa)
%                 temp = temperture (C)
%                 humi = relative humidity (%)
%            default:pre,temp=standard atmosphere,humi=50%
% [argout] : trop = tropospheric delay (m)
% [note]   : Reference: Montenbruck Satellite Orbits p.223, 1/É…'^2=0
% [version]: $Revision: 3 $ $Date: 06/07/08 14:22 $
%            Copyright(c) 2004-2006 by T.Takasu, all rights reserved
% [history]: 03/12/12   0.1  new
%-------------------------------------------------------------------------------
if nargin<4|isempty(met_prm) % standard atmosphere
    met_prm=[1013.25*(1-2.2557E-5*gpos(3))^5.2568,15.0-6.5E-3*gpos(3),50];
end
T=met_prm(2)+273.16; % temperture(K)
e=6.108*met_prm(3)/100*exp((17.15*T-4684)/(T-38.45)); % partial press. of water(hPa)
n=[77.624*met_prm(1)/T,371900*e/T^2-12.92*e/T];
h=11385*[met_prm(1)/n(1),(1255/T+0.05)*e/n(2)]; % height ot troposphere(m)
for j=1:2, dr(j)=TropPoly(azel(2),h(j)); end
trop=9.9842E-7*sum(n.*dr);

% tropospheric polynominal -----------------------------------------------------
function dr=TropPoly(el,h)
Re=6378137;
r=sqrt((Re+h)^2-(Re*cos(el))^2)-Re*sin(el);
a=-sin(el)/h; b=-cos(el)^2/(2*h*Re);
alph=[1,4*a,6*a*a+4*b,4*a*(a*a+3*b),a^4+12*a*a*b+6*b*b,4*a*b*(a*a+3*b),...
      b*b*(6*a*a+4*b),4*a*b^3,b^4]; % Éøij
dr=sum(alph.*(r.^(1:9))./(1:9));
