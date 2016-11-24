function [x,P,vz,stat]=filtersrcf(x,P,dz,G,R,nsig)
%-------------------------------------------------------------------------------
% [system] : GpsTools
% [module] : square root convergence filter
% [func]   : calculate measurement update rule by square root convergence filter
% [argin]  : x   = prefit states
%            P   = prefit covariences matrix
%            dz  = residuals
%            G   = jacobian matrix
%            R   = measurement noise corvariences matrix
%           (nsig) = outlier threshold (sigma,0:no exclusion)
% [argout] : xu  = postfit states
%            Pu  = postfit covarience matrix
%            vz  = valid residulas (1:valid,0:invalid/excluded)
%            stat= status (0:ok,<0:error,singular)
% [note]   :
% [version]: $Revision: 16 $ $Date: 2008-12-12 15:49:30 +0900 (金, 12 12 2008) $
% [history]: 04/07/22   0.1  new
%            08/11/21   0.2  mex-file with intel MKL -> m-file (gt_0.6.4)
%-------------------------------------------------------------------------------
if nargin<6, nsig=0; end
C=chol(P)';
V=C'*G';
S=V'*V+R;
if nsig>0, vz=dz.^2<=diag(S)*nsig^2; else vz=ones(length(dz),1); end
j=find(vz);
if isempty(j) % no data
    stat=-9999;
elseif cond(S(j,j))>1E14 % singular
    stat=-1;
else
    VS=V(:,j)/S(j,j);
    K=C*VS;
    x=x+K*dz(j);
    C=C*chol(eye(length(x))-VS*V(:,j)')';
    P=C*C';
    stat=0;
end
