function derp=erpvar(tutc,utc_tai)
%-------------------------------------------------------------------------------
% [system] : GpsTools
% [module] : earth rotation parameter daily and sub-daily variation
% [func]   : calculate daily and sub-daily variation of erp by earth tides
% [argin]  : tutc    = date/time (mjd-utc)
%            utc_tai = UTC-TAI (sec)
% [argout] : derp = variation of erp [dxp,dyp,dut]
%                dxp,dyp = xp,yp delta (rad)
%                dut     = UT1-UTC delta (sec)
% [note]   : reference : IERS Conventions 1996, Ch.8 Ray model
% [version]: $Revision: 12 $ $Date: 2008-11-25 10:02:15 +0900 (ÁÅ´, 25 11 2008) $
%            Copyright(c) 2004-2006 by T.Takasu, all rights reserved
% [history]: 04/11/19  0.1  new
%-------------------------------------------------------------------------------
mas2rad=4.8481368110953598E-9;
Nc=[
-1 0 -2 0 -2 1 -90 26.868  0.02  0.05 -0.026  0.006 -0.006 -0.026 % Q1
 0 0 -2 0 -2 1 -90 25.819  0.12  0.16 -0.133  0.049 -0.049 -0.133 % O1
 0 0 -2 2 -2 1 -90 24.066  0.03  0.05 -0.050  0.025 -0.025 -0.050 % P1
 0 0  0 0  0 1  90 23.935  0.09  0.18 -0.152  0.078 -0.078 -0.152 % K1
-1 0 -2 0 -2 2   0 12.658 -0.04 -0.02 -0.057 -0.013  0.011  0.033 % N2
 0 0 -2 0 -2 2   0 12.421 -0.16 -0.07 -0.330 -0.028  0.037  0.196 % M2
 0 0 -2 2 -2 2   0 12.000 -0.08  0.00 -0.145  0.064  0.059  0.087 % S2
 0 0  0 0  0 2   0 11.967 -0.02  0.00 -0.036  0.017  0.018  0.022 % K2
];
tt=tutc+(-utc_tai+32.184)/86400;
theta=(100.4606184+360.98564737*(tt-51544.5))*pi/180;
args=Nc(:,1:7)*[EphSunMoon(tt);theta;pi/180];
sina=sin(args);
cosa=cos(args);
dut=sum(Nc(:, 9).*sina+Nc(:,10).*cosa)*1E-4;    % UT1-UT1D(sec)
dxp=sum(Nc(:,11).*sina+Nc(:,12).*cosa)*mas2rad; % É¢x(rad)
dyp=sum(Nc(:,13).*sina+Nc(:,14).*cosa)*mas2rad; % É¢y(rad)
derp=[dxp,dyp,dut];

% sun/moon arguments -----------------------------------------------------------
function f=EphSunMoon(tt)
sec2rad=4.8481368110953598E-6;
Fc=[
134.96340251 1717915923.2178  31.8792  0.051635 -0.00024470 % l
357.52910918  129596581.0481  -0.5532  0.000136 -0.00001149 % l'
 93.27209062 1739527262.8478 -12.7512 -0.001037  0.00000417 % F
297.85019547 1602961601.2090  -6.3706  0.006593 -0.00003169 % D
125.04455501   -6962890.2665   7.4722  0.007702 -0.00005939 % OMG
];
t=(tt-51544.5)/36525;
f=mod(Fc*[pi/180;t*sec2rad;t^2*sec2rad;t^3*sec2rad;t^4*sec2rad],2*pi);
