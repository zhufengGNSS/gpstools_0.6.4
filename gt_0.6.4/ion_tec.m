function ion=ion_tec(td,time,azel,gpos,iondir,ionsrc)
%-------------------------------------------------------------------------------
% [system] : GpsTools
% [module] : ionospheric model - global tec map
% [func]   : calculate ionospheric delay by global tec map
% [argin]  : td   = date (mjd-gpst)
%            time = time vector (sec)
%            azel = satellite azimath/elevation angle (rad)
%                azel(n,:) = time(n) az/el [az,el]
%            gpos = station latilute/longitude/height(deg,m) [lat,lon,h]
%           (iondir) = tec data directory
%           (ionsrc) = tec data source
%                'igs'  : IGS tec map final
%                'igr'  : IGS tec map rapid
% [argout] : ion = ionospheric delay (L1) (m)
%                ion(n,1) = time(n) ion-delay
% [note]   :
% [version]: $Revision: 12 $ $Date: 2008-11-25 10:02:15 +0900 (火, 25 11 2008) $
%            Copyright(c) 2004-2006 by T.Takasu, all rights reserved
% [history]: 04/11/02  0.3  new
%-------------------------------------------------------------------------------
if nargin<5, iondir=''; end
if nargin<6, ionsrc=''; end
if isempty(ionsrc), ionsrc='igs'; end
ion=repmat(nan,length(time),1);
for n=1:length(time)
    dv=mjdtocal(td,time(n)); dirs=gfilepath(iondir,'',dv(1:3));
    ion(n)=IonDelay(td,time(n),azel(n,:),gpos,dirs,ionsrc);
end

% calculate ionospheric delay (L1) ---------------------------------------------
function ion=IonDelay(td,t,azel,gpos,iondir,ionsrc)
persistent tdd time tec lats lons hgts rb, if isempty(tdd), tdd=0; end
f1=1.57542E9; ion=nan; if azel(2)<0, return, end

tdr=floor(td+t/86400);
if tdr~=tdd
    [epoch,time,tec,lats,lons,hgts,rb]=ReadTecMap(tdr,iondir,ionsrc);
    if isempty(tec), return, end
    [tdd,ts]=caltomjd(epoch);
end
%[gpos,azel]=GeodToGeoc(gpos,azel);
[lat,lon,z]=IonPos(gpos,azel,hgts,rb);
vtec=InterpTec((td-tdd)*86400+t,lat,lon,tec,time,lats,lons);
ion=0; for n=1:length(vtec), ion=ion+40.3E16/cos(z(n))/f1^2*vtec(n); end

% read tec map ----------------------------------------------------------------
function [epoch,time,tec,lats,lons,hgts,rb]=ReadTecMap(td,iondir,ionsrc)
epoch=[]; time=[]; tec=[]; lats=[]; lons=[]; hgts=[]; rb=0;

switch ionsrc
case {'igs','igr'}
    dt=mjdtocal(td);
    day=caltomjd(dt(1:3))-caltomjd([dt(1),1,1])+1;
    file=sprintf('%sg%03d%1d.%02di',ionsrc,day,0,mod(dt(1),100));
    [epoch,time,tec,rms,lats,lons,hgts,rb]=readionex(fullfile(iondir,file));
    lats=lats*pi/180;
    lons=lons*pi/180;
otherwise,
end

% geodetic latitude->geocentric latitude  --------------------------------------
function [gpos,azel]=GeodToGeoc(gpos,azel)
[epos,Eg]=geodtoecef(gpos);
gpos(1)=asin(epos(3)/norm(epos))*180/pi;
[epos,Ec]=geodtoecef(gpos);
e=Ec*Eg'*[cos(azel(2))*sin(azel(1));cos(azel(2))*cos(azel(1));sin(azel(2))];
azel=[atan2(e(1),e(2)),asin(e(3))];

% pieas point position/zenith angle --------------------------------------------
function [lat,lon,z]=IonPos(gpos,azel,hgts,rb)
alpha=0.9782; % modified single-layer model
for n=1:length(hgts)
    latr=gpos(1)*pi/180; lonr=gpos(2)*pi/180; zr=pi/2-azel(2);
    z(n,1)=asin(rb/(rb+hgts(n))*sin(alpha*zr));
    g=zr-z(n,1);
    lat(n,1)=asin(cos(g)*sin(latr)+sin(g)*cos(latr)*cos(azel(1)));
    lon(n,1)=lonr+asin(sin(g)*sin(azel(1))/cos(lat(n,1)));
    if lon(n,1)<=-pi, lon(n,1)=lon(n,1)+2*pi;
    elseif lon(n,1)>pi, lon(n,1)=lon(n,1)-2*pi; end
end

% interporate tec map ----------------------------------------------------------
function vtec=InterpTec(t,lat,lon,tec,time,lats,lons)
vtec=repmat(nan,size(tec,3),1);
i=floor((length(time)-1)*(t-time(1))/(time(end)-time(1)))+1;
for n=1:size(tec,3)
    lon1=lon(n)+2*pi/86400*(t-time(i));
    lon2=lon(n)+2*pi/86400*(t-time(i+1));
    if ~isnan(lon1)&~isnan(lon2)
        e1=interp2(lons,lats,tec(:,:,n,i),lon1,lat(n),'*linear');
        e2=interp2(lons,lats,tec(:,:,n,i+1),lon2,lat(n),'*linear');
        vtec(n,1)=1/(time(i+1)-time(i))*((time(i+1)-t)*e1+(t-time(i))*e2);
    end
end
