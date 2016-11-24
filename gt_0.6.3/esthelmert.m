function tr=esthelmert(posr,poss,nsig)
%-------------------------------------------------------------------------------
% [system] : GpsTools
% [module] : estimate Helmert transformation parameters
% [func]   : estimate Helmert transformation parameters by square root adjustment
% [argin]  : posr = reference positions   [x1,y1,z1;x2,y2,z2;..] (m)
%            poss = transformed positions [x1,y1,z1;x2,y2,z2;..] (m)
%           (nsig)= outlier threshold (sigma,0:no outler exclusion) (defalut:0)
% [argout] : tr  = helmert transformation parameters(error:[0;0;0;0;0;0;0]
%                tr(1) = dx (m)
%                tr(2) = dy (m)
%                tr(3) = dz (m)
%                tr(4) = Rx (mas)
%                tr(5) = Ry (mas)
%                tr(6) = Rz (mas)
%                tr(7) = scale (ppb)
% [note]   :
% [version]: $Revision: 4 $ $Date: 06/07/08 14:11 $
%            Copyright(c) 2004-2006 by T.Takasu, all rights reserved
% [history]: 05/03/05   0.1  new
%            05/06/27   0.2  add outlier exclusion
%-------------------------------------------------------------------------------
if nargin<3, nsig=0; end
i=find(~isnan(posr(:,1))&~isnan(poss(:,1)));
while 1
    [tr,res]=EstTrMat(posr(i,:),poss(i,:));
    resp=sum(res.^2,2);
    [maxr,j]=max(resp);
    if ~any(tr)|nsig==0|maxr<=mean(resp)*nsig^2, break, else i(j)=[]; end
end

% estimate transformation matrix -----------------------------------------------
function [tr,res]=EstTrMat(posr,poss)
m2r=pi/180/3600*1E-3; % mas->rad
AA=zeros(7); Ay=zeros(7,1);
for n=1:size(posr,1)
    r=posr(n,:)';
    A=[eye(3),[0,r(3),-r(2);-r(3),0,r(1);r(2),-r(1),0]*m2r,r*1E-9];
    AA=AA+A'*A;
    Ay=Ay+A'*(poss(n,:)'-r);
end
if cond(AA)<=1E15, tr=AA\Ay;
else disp('warning : esthelmert error'); tr=zeros(7,1); end
res=poss-helmert(posr,tr);
