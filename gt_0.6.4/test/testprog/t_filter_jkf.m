function t_filter_jkf
%-------------------------------------------------------------------------------
% カルマンフィルタ(Joseph)試験 2004/01/08 by TTAKA
%-------------------------------------------------------------------------------
nx=1000; nz=500;

dz=rand(nz,1);
x=rand(nx,1);
P=rand(nx); P=P*P';
G=full(sprand(nz,nx,0.001));
sig=rand(nz,1);

[xe,Pe]=filter_jkf(dz,x,P,G,sig);
[xr,Pr]=filter_jkfr(dz,x,P,G,sig);
max(max(xe-xr))
max(max(Pe-Pr))

% カルマン(Joseph)
function [xu,Pu]=filter_jkfr(dz,x,P,G,sig)
G=sparse(G);
R=diag(sig.^2);
K=P*G'*inv(G*P*G'+R);
xu=x+K*dz;
D=eye(length(x))-K*G;
Pu=D*P*D'+K*R*K';
