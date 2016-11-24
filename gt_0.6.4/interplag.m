function [yi,dyi]=interplag(x,y,xi,nmax)
%-------------------------------------------------------------------------------
% [system] : GpsTools
% [module] : lagrange interpolation
% [func]   : lagrange interpolation by sliding segment
% [argin]  : x    = x
%            y    = f(x)
%            xi   = interpolated x values
%           (nmax)= degree of interpolation polynominal (default:8)
% [argout] : yi   = f(xi)
%            dyi  = f'(xi)
% [note]   :
% [version]: $Revision: 12 $ $Date: 2008-11-25 10:02:15 +0900 (火, 25 11 2008) $
%            Copyright(c) 2004-2006 by T.Takasu, all rights reserved
% [history]: 05/09/20   0.1  new
%            06/03/31   0.2  improve performance
%-------------------------------------------------------------------------------
if nargin<4, nmax=8; end
if nmax>=length(x), nmax=length(x)-1; end
if nmax<1, error('nmax too small'); end
yi=repmat(nan,length(xi),size(y,2)); dyi=yi;
ns=floor((nmax+1)/2);

% interpolate head segment
i=find(xi<x(ns));
if ~isempty(i)
    [yi(i,:),dyi(i,:)]=interpseg(x(1:nmax+1),y(1:nmax+1,:),xi(i));
end
% interpolate intermediate segments
for n=ns:length(x)-nmax+ns-1
    i=find(x(n)<=xi&xi<x(n+1));
    if ~isempty(i)
        [yi(i,:),dyi(i,:)]=interpseg(x((1:nmax+1)+n-ns),y((1:nmax+1)+n-ns,:),xi(i));
    end
end
% interpolate tail segments
i=find(x(end-nmax+ns)<=xi);
if ~isempty(i)
    [yi(i,:),dyi(i,:)]=interpseg(x(end-nmax:end),y(end-nmax:end,:),xi(i));
end

% interpolate segment ----------------------------------------------------------
function [yi,dyi]=interpseg(x,y,xi)
nx=length(x); ni=length(xi);
xi=xi(:); yi=zeros(ni,size(y,2)); zi=zeros(ni,size(y,2)); dx=1E-3;
for n=1:nx
    xx=x; xx(n)=[]; a=prod(x(n)-xx);
    b=xi*ones(1,nx-1)-ones(ni,1)*xx;
    yi=yi+prod(b,2)*y(n,:)/a;
    zi=zi+prod(b+dx,2)*y(n,:)/a;
end
dyi=(zi-yi)/dx; % derivative
