function T_GenOrbit
%-------------------------------------------------------------------------------
% ãOìπê∂ê¨ï]âø
%-------------------------------------------------------------------------------
global eoptype, eoptype=0;

tm1=CalToMjd([2008,1,1,0,0,0]);
[U,P,N,gmst]=Ecef(tm1);

%tm0=mjd([2008,1,1,0,0,0]);
%P,Precession(tm0)'
%N,Nutation(tm0)'
%U,TrsToCrs(tm0)'
%gmst,GMST(tm0)

load('QZS01_20080101_EPH')
PlotOrbit(toe,time,ephem,'.')
%VerifyOrbit(toe,time,ephem)

% ãOìπÉvÉçÉbÉg -----------------------------------------------------------------
function PlotOrbit(toe,time,ephem,style)
global eoptype, eoptype=0;
t0=Mjd(toe);
r0=LLHtoXYZ(36*pi/180,140.1*pi/180,0);
for n=1:length(time)
   rs=TRStoCRS(t0+time(n)/86400)'*ephem(n,1:3)';
   [lat(n),lon(n)]=XYZtoLLH(rs);
   [az(n),el(n)]=XYZtoAZEL(rs,r0);
end
%WorldMap(lat'*180/pi,lon'*180/pi,style);
SkyMap(az*180/pi,el*180/pi,'b.');

%subplot(2,1,1), plot(time/3600,el)
%xlabel('hour'), ylabel('elevation(deg)')
%axis([0 time(end)/3600 0 90])
%xlabel('hour'), ylabel('elevation(deg)')
%set(gca,'xtick',0:24,'ytick',0:15:90), grid on
%subplot(2,1,2), plot(time/3600,az)
%axis([0 time(end)/3600 0 360])
%xlabel('hour'), ylabel('azimuth(deg)')
%set(gca,'xtick',0:24,'ytick',0:90:360), grid on

% ãOìπî‰är ---------------------------------------------------------------------
function VerifyOrbit(toe,time,ephem)
global Pflag geonmax SrpPar SatMass eoptype

Pflag=[1,1,1,1,0,0,0,0]; geonmax=8; SrpPar=[1.3 18.4]; SatMass=590; eoptype=0;
toe=[2008,1,1,0,0,0];
ele=[42164000.000,0.100000000000,45.0000000,326.0294000,270.0000000,0.0000000];
ele(3:6)=ele(3:6)*pi/180;
[r0,v0]=ELEtoXYZ(ele(1),ele(2),ele(3),ele(4),ele(5),ele(6));
[t,r,v]=PredOrbit(Mjd(toe),time,r0,v0,0);

for n=1:length(time)
    dr(n)=norm(r(n,:)-ephem(n,1:3));
    dv(n)=norm(v(n,:)-ephem(n,4:6));
end
plot(time,dr)
axis([0 3600 0 0.01]), grid on
xlabel('sec'), ylabel('pos error(m)')
