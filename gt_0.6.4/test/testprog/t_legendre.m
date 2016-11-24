function t_legendre

[p0,dp0]=CalcLegendre(6,0.3);
[p1,dp1]=LegendreFunc(0.3,6);
p1-p0
dp1-dp0

% ³‹K‰»ƒ‹ƒWƒƒƒ“ƒhƒ‹”†ŠÖ” -----------------------------------------------------
% [in]  : x  = x
%         nmax = Å‘åŽŸ”
% [out] : pnm = Pnm(x) (pnm(n,m+1)=Pnm)
%         dpnm = Pnm'(x) (dpnm(n,m+1)=Pnm')
%-------------------------------------------------------------------------------
function [pnm,dpnm]=LegendreFunc(x,nmax)
persistent npnm, if isempty(npnm), npnm=NormalLegendreCoef; end
pnm=zeros(nmax+1,nmax+1); dpnm=pnm;
pnm(1:2,1:2)=[1,0;x,sqrt(1-x*x)];
for n=3:nmax+1
    pnm(n,n)=(2*n-3)*sqrt(1-x*x)*pnm(n-1,n-1);
    dpnm(n,n)=(n-1)*x*pnm(n,n);
    for m=1:n-1
        pnm(n,m)=(x*(2*n-3)*pnm(n-1,m)-(n+m-3)*pnm(n-2,m))/(n-m);
        dpnm(n,m)=(n-1)*x*pnm(n,m)-(n+m-2)*pnm(n-1,m);
    end
end
pnm=npnm(1:nmax,1:nmax+1).*pnm(2:end,:);
dpnm=npnm(1:nmax,1:nmax+1).*dpnm(2:end,:)./(x*x-1);

% ³‹K‰»ŒW” -------------------------------------------------------------------
% [in]  : ‚È‚µ
% [out] : npnm = Pnm³‹K‰»ŒW”
%-------------------------------------------------------------------------------
function npnm=NormalLegendreCoef
npnm=zeros(10,10+1);
for n=1:10
    npnm(n,1)=sqrt(2*n+1);
    for m=1:n, npnm(n,m+1)=sqrt(factorial(n-m)*(4*n+2)/factorial(n+m)); end
end


function [P, DP] = CalcLegendre(nmax, x)
%-------------------------------------------------------------------------------
% CALCLEGENDRE(nmax, x)
%y ‹@”\ z: ³‹K‰»ƒ‹ƒWƒƒƒ“ƒhƒ‹”†ŠÖ”/”†“±ŠÖ”(Pnm(x),dPnm(x)/dx)ŒvŽZ
%y ˆø” z: nmax = Å‘åŽŸ”
%            x = ŠÖ”ˆø” x
%y–ß‚è’lz: P = Pnm(x)s—ñ(nmax x nmax+1) (P(n,m+1) = Pnm)
%            DP = P'nm(x)s—ñ(nmax x nmax+1) (DP(n,m+1) = P'nm)(È—ª‰Â)
%y “Á‹L z: DP‚Ì1ŽŸ‚ÍŒvŽZ‚µ‚È‚¢
%-------------------------------------------------------------------------------
persistent N
if isempty(N), N = CalcNormCoef(12); end         % N=[Nmn] ³‹K‰»ŒW”
dflg = nargout >= 2;

P = zeros(nmax+1,nmax+1);
if dflg, DP = zeros(nmax+1,nmax+1); end
xx1 = 1 - x .* x;
P(1:2,1) = [1; x];
for n = 2:nmax
    % [P10(x) P20(x) P30(x) ... Pn0(x)]ŒvŽZ
    P(n+1,1) = (x .* (2 .* n - 1) .* P(n,1) - (n - 1) .* P(n-1,1)) ./ n;
    
    % ‘Q‰»Ž®‚Å[P20'(x) P30'(x) ... Pn0'(x)]*(x^2-1)ŒvŽZ
    if dflg, DP(n+1,1) = n .* (x .* P(n+1,1) - P(n,1)); end
end
for m1 = 2:nmax+1 % m1=m+1
    % [P11(x) P22(x) ... Pnn(x)], [P11'(x) P22'(x) ... Pnn'(x)]*(x^2-1)ŒvŽZ
    P(m1,m1) = prod(2.*m1-3:-2:1) .* xx1 .^ ((m1 - 1) ./ 2);
    if dflg, DP(m1,m1) = (m1 - 1) .* x .* P(m1,m1); end
    
    % [P21(x) P31(x) ... P32(x) P42(x) ... P43(x) P53(x) ...]
    % [P21'(x) P31'(x) ... P32'(x) P42'(x) ... P43'(x) P53'(x) ...]*(x^2-1)ŒvŽZ
    for n = m1:nmax
        P(n+1,m1) = ...
        (x .* (2 .* n - 1) .* P(n,m1) - (n + m1 - 2) .* P(n-1,m1)) ./ (n - m1 + 1);
        if dflg, DP(n+1,m1) = n .* x .* P(n+1,m1) - (n + m1 - 1) .* P(n,m1); end
    end
end
% ƒ‹ƒWƒƒƒ“ƒhƒ‹”†ŠÖ”,”†“±ŠÖ”³‹K‰»(->P,DP)
P = N(1:nmax,1:nmax+1) .* P(2:end,:);
if dflg, DP = - N(1:nmax,1:nmax+1) .* DP(2:end,:) ./ xx1; end


function N = CalcNormCoef(nmax)
%-------------------------------------------------------------------------------
% CALCNORMCOEF(nmax)
%y ‹@”\ z: ƒ‹ƒWƒƒƒ“ƒhƒ‹”†ŠÖ”³‹K‰»ŒW”ŒvŽZ
%y ˆø” z: nmax = Å‘åŽŸ”
%y–ß‚è’lz: N = Nnms—ñ(nmax x nmax+1) (N(n,m+1)=Nnm)
%y “Á‹L z:
%-------------------------------------------------------------------------------
N = zeros(nmax, nmax+1);
for n = 1:nmax
    for m = 0:n
        % ³‹K‰»ŒW”NnmŒvŽZ(->N)
        if m <= 0, N(n,1) = sqrt(2 .* n + 1);
        else
            N(n,m+1) = sqrt(factorial(n - m) .* (4 .* n + 2) ./ factorial(n + m));
        end
    end
end

