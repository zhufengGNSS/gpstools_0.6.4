function [data,gprm,layer]=readgpv(td,ts,type,gpvdir,gpvsrc,layer,ft)
%-------------------------------------------------------------------------------
% [system] : GpsTools
% [module] : read gpv data
% [func]   : read gpv data
% [argin]  : td,ts   = date (mjd),time (sec)
%            type    = parameter type
%                'pmsl' : mean sea level pressure (hPa)
%                'geoh' : geopotential height (m)
%                'temp' : temperture (C)
%                'wind' : wind (m/s)
%                'vvel' : vertical velocity (hPa/s)
%                'humi' : relative humidity (%)
%                'rain' : rain (mm/s)
%                'clou' : cloud (rate)
%                'geoid': geoid height (m)
%           (gpvdir) = gpv data directory (default:current)
%           (gpvsrc) = gpv data source (default:'rsm')
%                    'rsm' : JMA RSM model
%                    'msm' : JMA MSM model
%                    'gsm' : JMA GSM model
%                    'rso' : JMA RSM model online
%                    'mso' : JMA MSM model online
%                    'gso' : JMA GSM model online
%                    'egm' : EGM96 Geoid
%                    'gsi' : GSI Geoid
%           (layer)  = parameter layer (hPa) (nan;all,0:surface)
%           (ft)     = forecast time (hr)
% [argout] : data  = gpv data
%                data(:,:,n,m) : layer(n) gpv data (m=1:winu,m=2:winv)
%            gprm  = grid parameters
%              gprm.nv   : number of vertical coordinate parameters
%              gprm.pv   : location of vertical coordinate parameters
%              gprm.type : data representation type
%                            0 : latitude/longitude grid (eq-cylindrical)
%                            1 : mercator projection
%                            3 : lambert conformal projection
%                            5 : polar stereographic projection
%                           13 : oblique lambert conformal projection
%                           50 : spherical harmonic coefficients
%                           90 : space view of orthgraphic grid
%              gprm.nx,ny : no. of points along lat/lon or x/y-axis
%              gprm.rcflg : resolution and component flag
%                           rcflg(1) : direction increment given
%                           rcflg(2) : earth radius
%                               0 : sphere re = 6367.47km
%                               1 : oblate spheriod re=6378.16km f=1/297.0
%                           rcflg(3) : uv component of vector
%                               0 : relative to east/north
%                               1 : relative to defined grid x/y
%              gprm.smode : scanning mode
%                           0 : points scan in + direction
%                           1 : points scan in - direction
%              gprm.lat1,lon1 : lat/lon of first grid (deg or m)
%              gprm.lat2,lon2 : lat/lon of last grid (deg) (type=0)
%              gprm.dx,dy : lat/lon or x/y increment (deg or m)
%              gprm.lov   : orientation of grid (y-parallel lon.) (deg) (type!=0)
%              gprm.lati1,lati2 : lats which secant cone cuts the sphere (type!=0)
%              gprm.latp,lonp : lat/lon of southern pole (type!=0)
%              gprm.x0,y0 : projection base grid point
%            layer = parameter layers (hPa) (-1:surface)
% [note]   :
% [version]: $Revision: 8 $ $Date: 06/07/18 15:31 $
%            Copyright(c) 2004-2006 by T.Takasu, all rights reserved
% [history]: 05/05/02  0.1  new
%-------------------------------------------------------------------------------
if nargin<4, gpvdir=''; end
if nargin<5, gpvsrc=''; end
if nargin<6, layer=nan; end
if nargin<7, ft=0; end
if isempty(gpvsrc), gpvsrc='rsm'; end
[data,gprm,layer]=ReadGpvData(td,ts,type,gpvdir,gpvsrc,layer,ft);

% read gpv data --------------------------------------------------------------
function [data,gprm,layer]=ReadGpvData(td,ts,type,gpvdir,gpvsrc,layer,ft)
persistent file gdata, if isempty(file), file=''; end
data=[]; gprm=[];
switch gpvsrc
case {'gsm','rsm','msm','gso'} % JMA GSM/RSM/MSM OFFLINE, GSM ONLINE
    switch gpvsrc
    case 'gsm', name='GANAL'; ss={''};
    case 'rsm', name='RANAL'; ss={''};
    case 'msm', name='MANAL'; ss={''};
    case 'gso', name='GSM'; ss={'X024','X048'};
    end
    dt=mjdtocal(td,ts);
    if strcmp(gpvsrc,'gso'), fn=sprintf('%s%02d',name,dt(4));
    else fn=sprintf('%s_%04d%02d%02d%02d.grb',name,dt(1:4)); end
    f=gfilepath(gpvdir,fn,dt);
    if ~strcmp(f,file)
        gdata=[]; for n=1:length(ss), gdata=[gdata;readgrib([f,ss{n}])]; end
        if isempty(gdata), disp(['warning : grib data read error :',f]), return, end
        file=f;
    end
    if strcmp(type,'wind')
        [winu,gprmu,layeru]=GetData(gdata,33,layer,ft,strcmp(gpvsrc,'gso')); % win-u
        [winv,gprmv,layerv]=GetData(gdata,34,layer,ft,strcmp(gpvsrc,'gso')); % win-v
        if ~isempty(winu)&~isempty(winv)&gprmu.nx==gprmv.nx&gprmu.ny==gprmv.ny
            [layer,i,j]=intersect(layeru,layerv);
            data=cat(4,winu(:,:,i),winv(:,:,j));
            gprm=gprmu;
        end
    else
        types={'pmsl','geoh','temp','vvel','humi'};
        tts=[2,7,11,39,52];
        i=find(strcmp(type,types));
        if ~isempty(i)
            [data,gprm,layer]=GetData(gdata,tts(i),layer,ft,strcmp(gpvsrc,'gso'));
            if any(strcmp(type,{'pmsl','vvel'})) % Pa->hPa
                data=single(double(data)/100);
            elseif strcmp(type,'temp') % K->C
                data=single(double(data)-273.15);
            end
        end
    end
    if ~isempty(gprm)
        switch gpvsrc
        case 'rsm', gprm.x0=200; gprm.y0=185; % projection base grid
        case 'msm', gprm.x0=245; gprm.y0=205;
        otherwise, gprm.x0=0; gprm.y0=0;
        end
    end
case {'rso','mso'} % JMA RSM/MSM ONLINE
    switch gpvsrc
    case 'rso', name='RSM'; ss={'SFC024','SFC051','PLM024','PLM051','PMH024','PMH051'};
    case 'mso', name='MSM'; ss={'SFC018_1','SFC018_2','PLM018_1','PLM018_2','PMH018'};
    end
    dt=mjdtocal(td,ts);
    f=gfilepath(gpvdir,sprintf('%s%02d',name,dt(4)),dt);
    if ~strcmp(f,file)
        gdata=[]; for n=1:length(ss), gdata=[gdata;readdgrb([f,ss{n}])]; end
        if isempty(gdata), return, end
        file=f;
    end
    if strcmp(type,'wind')
        [winu,gprmu,layeru]=GetData2(gdata,23,layer,ft,gpvsrc); % wind u-component
        [winv,gprmv,layerv]=GetData2(gdata,24,layer,ft,gpvsrc); % wind v-component
        if ~isempty(winu)&~isempty(winv)&gprmu.nx==gprmv.nx&gprmu.ny==gprmv.ny
            [layer,i,j]=intersect(layeru,layerv);
            data=cat(4,winu(:,:,i),winv(:,:,j));
            gprm=gprmu;
        else
            layer=[];
        end
    else
        types={'pmsl','geoh','temp','vvel','humi','rain','clou'};
        tts=[1,102,4,42,13,49,225];
        i=find(strcmp(type,types));
        if ~isempty(i), [data,gprm,layer]=GetData2(gdata,tts(i),layer,ft,gpvsrc); end
    end
case {'gsi','egm'} % GSI Geoid2000
    if strcmp(gpvsrc,'gsi'), fn='geoid_gsi2000'; else fn='geoid_egm96'; end
    [dirs,f]=fileparts(which(mfilename));
    file=fullfile(dirs,'data',[fn,'.mat']);
    if exist(file), load(file); end

otherwise
    disp(['warning : no support gpv type : ',gpvsrc])
end

% get type data(grib) ----------------------------------------------------------
function [data,gprm,layer]=GetData(gdata,type,layer,ft,opt)
data=[]; gprm=[]; slayer=[];
if isempty(gdata), layer={}; return, end
pprm=[gdata.pprm]; lat=[-90,90]; lon=[-180,180];
for i=find([pprm.param]==type&[pprm.p1]==ft)
    if pprm(i).level(1)==1, level=0; else level=pprm(i).level(2); end
    if isnan(layer)|level==layer
        if opt % gsm global
            j=find(slayer==level);
            if isempty(j)
                gprm=gdata(i).gprm;
                gprm.nx=360/gprm.dx; gprm.ny=180/gprm.dy+1;
                gprm.lat1=90; gprm.lat2=-90; gprm.lon1=-180; gprm.lon2=180-gprm.dx;
                data=cat(3,data,repmat(nan,gprm.ny,gprm.nx));
                slayer=[slayer;level]; j=length(slayer);
            end
            n=-gdata(i).gprm.lat1/gprm.dy+(1:gdata(i).gprm.ny);
            m=(gdata(i).gprm.lon1+180)/gprm.dx+(1:gdata(i).gprm.nx);
            k=find(m>gprm.nx); m(k)=m(k)-gprm.nx;
            data(n,m,j)=flipud(gdata(i).data);
        elseif isempty(gprm)|(gprm.nx==gdata(i).gprm.nx&gprm.ny==gdata(i).gprm.ny)
            data=cat(3,data,flipud(gdata(i).data));
            gprm=gdata(i).gprm;
            slayer=[slayer;level];
        else
            disp('warning : not compatible data set');
        end
    end
end
if ~isempty(data), [slayer,i]=sortrows(slayer); data=data(:,:,i); end
layer=slayer;

% get type data(dgrb) ----------------------------------------------------------
function [data,gprm,layer]=GetData2(gdata,type,layer,ft,src)
data=[]; gprm=[]; slayer=[];
if isempty(gdata), layer=[]; return, end
for i=find([gdata.id]==type)
    j=find(gdata(i).ft==ft);
    if ~isempty(j)
        if isnan(layer)|gdata(i).press==layer
            if isnan(layer)&gdata(i).press==0 % surface
                data=cat(3,data,gdata(i).data(1:2:end,1:2:end,j)');
            else
                data=cat(3,data,gdata(i).data(:,:,j)');
            end
            slayer=[slayer;gdata(i).press];
        end
    end
end
if isempty(data), layer=[]; return, end
gprm.type=0; gprm.nx=size(data,2); gprm.ny=size(data,1);
gprm.rcflg=[0,0,0]; gprm.smode=[0,1];
gprm.lov=0; gprm.lati1=0; gprm.lati2=0; gprm.latp=0; gprm.lonp=0;
gprm.x0=0; gprm.y0=0;
switch src
case 'gso',
case 'rso'
    if layer==0, gprm.dx=0.25; gprm.dy=0.2;
    else         gprm.dx=0.5;  gprm.dy=0.4; end
    gprm.lon1=120+gprm.dx/2; gprm.lon2=gprm.lon1+gprm.dx*(gprm.nx-1);
    gprm.lat1= 50-gprm.dy/2; gprm.lat2=gprm.lat1-gprm.dy*(gprm.ny-1);
case 'mso'
    if layer==0, gprm.dx=0.125; gprm.dy=0.1;
    else         gprm.dx=0.25;  gprm.dy=0.2; end
    gprm.lon1=120 +gprm.dx/2; gprm.lon2=gprm.lon1+gprm.dx*(gprm.nx-1);
    gprm.lat1=47.6-gprm.dy/2; gprm.lat2=gprm.lat1-gprm.dy*(gprm.ny-1);
end
[layer,i]=sortrows(slayer); data=data(:,:,i);
