function pos=helmert(pos,tr)
%-------------------------------------------------------------------------------
% [system] : GpsTools
% [module] : Helmert transformation
% [func]   : Helmert transformation
% [argin]  : poss = position [x,y,z;...] (m)
%            tr   = helmert transformation parameters
%                tr(1) = dx (m)
%                tr(2) = dy (m)
%                tr(3) = dz (m)
%                tr(4) = Rx (mas)
%                tr(5) = Ry (mas)
%                tr(6) = Rz (mas)
%                tr(7) = delta s (scale) (ppb)
% [argout] : pos = transformed position [x',y',z';...] (m)
% [note]   :
% [version]: $Revision: 3 $ $Date: 06/07/08 14:12 $
%            Copyright(c) 2004-2006 by T.Takasu, all rights reserved
% [history]: 05/03/05   0.1  new
%            05/03/24   0.2  support vecters
%-------------------------------------------------------------------------------
t=eye(3)+[0,-tr(6),tr(5);tr(6),0,-tr(4);-tr(5),tr(4),0]'*pi/180/3600*1E-3;
for n=1:size(pos,1), pos(n,:)=pos(n,:)*t; end
pos=(1+tr(7)*1E-9)*pos+repmat(tr(1:3)',size(pos,1),1);
