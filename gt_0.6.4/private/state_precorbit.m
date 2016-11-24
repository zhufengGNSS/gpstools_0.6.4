function [x,phi]=state_precorbit(t0,t,x0,sat,oprm)
%-------------------------------------------------------------------------------
% [system] : GpsTools
% [module] : state transition model - precise satellite orbit model
% [func]   : propagate states tempolary by precise satellite orbit model
% [argin]  : t0    = initial time (mjd-utc)
%            t     = time relative to t0 (sec)
%            x0    = initial satellite position/velocity [r0;v0;p0;srpprm] (eci)
%           (sat)  = satellite
%           (oprm) = satellite orbtit parameters
%                    (default:8,'grv_sunmoon','srp_code','','','','',30,'iers1996')
%                oprm.g_nmax    = degree of geopotential
%                oprm.p_plgrv   = solar/planetary potentials
%                oprm.p_solarpr = solar radiation pressure model
%                oprm.p_eclipse = eclipse model
%                oprm.p_tidal   = earth tides
%                oprm.p_relativ = relativistic effects
%                oprm.p_deltav  = thrust forces
%                oprm.tstep     = integration step time (sec)
%                oprm.satprm    = satellite parameters
%                oprm.nutmodel  = nutation model
% [argout] : x     = satellite position/velocity at t [r;v;p;srpprm] (eci)
%           (phi)  = state transition/sensitibity matrix phi(t,t0)=dx/dx0
% [note]   : dprm/dr=0,dprm/dv=0
% [version]: $Revision: 12 $ $Date: 2008-11-25 10:02:15 +0900 (火, 25 11 2008) $
%            Copyright(c) 2004-2006 by T.Takasu, all rights reserved
% [history]: 04/11/19  0.1  new
%            05/03/23  0.2  add oprm.satprm
%            05/06/13  0.3  add argin t of prm_gpssats()
%                           add argin oprm.nutmodel
%                           change default parameters
%                           change dx value in differential approx.
%-------------------------------------------------------------------------------
global ephpdir
if isempty(ephpdir)
    [dirs,f]=fileparts(which(mfilename));
    ephpdir=fullfile(dirs,'..','data');
end
if nargin<4, sat=''; end
if nargin<5, oprm=[]; end
if isempty(oprm)
    oprm.g_nmax   =8;
    oprm.p_plgrv  ='grv_sunmoon';
    oprm.p_solarpr='srp_code';
    oprm.p_eclipse='';
    oprm.p_tidal  ='';
    oprm.p_relativ='';
    oprm.p_deltav ='';
    oprm.tstep    =30;
    oprm.satprm   =[0,975,12,-90,-0.7,0.4,1.0]; % [type,mass,x-sec,D0,Y0,B0,Z0]
    oprm.nutmodel ='iers1996';
end
if length(x0)>6, sprm=x0(7:end)'; else sprm=DefaultSrp(sat,oprm.p_solarpr,t0); end

% set satellite orbit parameters
global satobt_toe satobt_prm satobt_srp satobt_ecl satobt_satprm satobt_tstep
satobt_toe=t0;
satobt_prm={oprm.g_nmax,oprm.p_plgrv,oprm.p_solarpr,oprm.p_eclipse,...
            oprm.p_tidal,oprm.p_relativ,oprm.p_deltav,oprm.nutmodel};
satobt_srp=sprm;
satobt_ecl=0;
satobt_tstep=oprm.tstep;
satobt_satprm=oprm.satprm;

% solve ordinary differential eq. of satellite orbit
x=x0; x(1:6)=satorbit_s(t,x0(1:6));

% partial derivatives
if nargout>=2
    phi=eye(length(x0));
    phi(1:6,1:6)=DxDx0(t0,t,x0,x);
    if length(x0)>6
        phi(1:6,7:end)=DxDp(t0,t,x0,x,oprm.p_solarpr);
    end
end

% state transition matrix ------------------------------------------------------
function dxdx0=DxDx0(t0,t,x0,x)
dx=[1,1,1,1E-4,1E-4,1E-4];
dxdx0=zeros(6,6);
for n=1:6 % differential approx.
    xx0=x0(1:6); xx0(n)=x0(n)+dx(n);
    dxdx0(:,n)=(satorbit_s(t,xx0)-x(1:6))/dx(n);
end

% sensitivity matrix with srp parameters ---------------------------------------
function dxdp=DxDp(t0,t,x0,x,model)
global satobt_srp
switch model
case 'srp_simple', dp=[1E-5];                    % [Cr*A/mass]
case 'srp_rock4',  dp=[1E-3,1E-3];               % [scale,ybias]
case 'srp_gspm',   dp=[1E-3,1E-3];               % [scale,ybias]
case 'srp_gspmm',  dp=[1E-3,1E-3,1E-3,1E-3];     % [scale,ybias,along,cross]
case 'srp_code',   dp=[1,1,1];                   % [D0,Y0,B0]
case 'srp_code2',  dp=[1,1,1,1];                 % [D0,Y0,B0,Z0]
case 'srp_code3',  dp=[1,1,1,1,1,1];             % [D0,Y0,B0,Z0,X10,X30]
end
dxdp=zeros(6,length(x)-6);
for n=1:length(dp) % differential approx.
    satobt_srp=x0(7:end)'; satobt_srp(n)=satobt_srp(n)+dp(n);
    dxdp(:,n)=(satorbit_s(t,x0(1:6))-x(1:6))/dp(n);
end

% get default srp parameters ---------------------------------------------------
function sprm=DefaultSrp(sat,model,t)
satprm=prm_gpssats(t);
i=find(strcmp(satprm(:,1),sat));
if isempty(i), srp=[0,975,12,-91.0,-0.9,1.3,0.9]; else srp=[satprm{i,2:end}]; end
switch model
case 'srp_simple', sprm=[1.3*srp(3)/srp(2)]; % [Cr*A/mass]
case 'srp_rock4',  sprm=[1,0];               % [scale,ybias]
case 'srp_gspm',   sprm=[1,0];               % [scale,ybias]
case 'srp_gspmm',  sprm=[1,0,0,0];           % [scale,ybias,along,cross]
case 'srp_code',   sprm=[srp(4:6)];          % [D0,Y0,B0]
case 'srp_code2',  sprm=[srp(4:7)];          % [D0,Y0,B0,Z0]
case 'srp_code3',  sprm=[srp(4:7),0,0];      % [D0,Y0,B0,Z0,X10,X30]
otherwise sprm=[];
end
