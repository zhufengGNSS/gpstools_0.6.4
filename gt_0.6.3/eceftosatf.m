function E=eceftosatf(state)
%-------------------------------------------------------------------------------
% [system] : GpsTools
% [module] : transform ecef to satellite orbit coordinate
% [func]   : transform ecef to satellite orbit coordinate
% [argin]  : state = satellite position/velocity[x;y;z;vx;vy;vz](m) (ecef)
% [argout] : E = transformation matrix ecef to satellite orbit coordinate
%                [r_radial;r_along-track;r_cross-track]=E*r
% [note]   :
% [version]: $Revision: 3 $ $Date: 06/07/08 14:10 $
%            Copyright(c) 2004-2006 by T.Takasu, all rights reserved
% [history]: 05/03/06   0.1  new
%-------------------------------------------------------------------------------
omge=7.2921151467E-5;
pos(:,1)=state(1:3);
vel(:,1)=state(4:6)+Cross([0;0;omge],pos);
crt=Cross(pos,vel);
alt=Cross(crt,pos);
E=[pos/norm(pos),alt/norm(alt),crt/norm(crt)]';

% cross (3x1) ------------------------------------------------------------------
function z=Cross(x,y)
z=[x(2)*y(3)-x(3)*y(2);x(3)*y(1)-x(1)*y(3);x(1)*y(2)-x(2)*y(1)];
