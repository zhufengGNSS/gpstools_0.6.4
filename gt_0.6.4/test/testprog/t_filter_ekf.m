function t_filter_ekf
%-------------------------------------------------------------------------------
% カルマンフィルタ(標準)試験 2004/01/08 by TTAKA
%-------------------------------------------------------------------------------
nx=633; nz=499;
%nx=72; nz=44;

dz=rand(nz,1);
x=rand(nx,1);
P=rand(nx); P=P*P';
G=full(sprand(nz,nx,0.001));
sig=rand(nz,1);

[xe,Pe]=filter_ekf(x,P,G,dz,sig);
[xr,Pr]=filter_ekfr(x,P,G,dz,sig);
max(max(xe-xr))
max(max(Pe-Pr))

% カルマン更新則
function [xu,Pu]=filter_ekfr(x,P,G,dz,sig)
G=sparse(G);
F=P*G';
K=F*inv(G*F+diag(sig.^2));
xu=x+K*dz;
Pu=P-K*G*P;
