function gencoast(file,coast,dirs,exlon,exlat,minland,minlake,delriver,mindp)
%-------------------------------------------------------------------------------
% [system] : GpsTools
% [module] : generate coast line data
% [func]   : read coast line database and generate coast line mat-file
% [argin]  : file    = output file path
%            coast   = input coast line database
%           (dirs)   = input data directory
%           (exlon)  = extent of longitude
%           (exlat)  = extent of latitude
%           (minland)= minimum area of land (km^2)
%           (minlake)= minimum area of lake (km^2)
%           (delriver)= delete river
%           (mindp)  = reduce polygon points
% [argout] : none
% [note]   : gshhs db
% [version]: $Revision: 1 $ $Date: 06/07/08 8:57 $
% [history]: 04/10/29  0.1  new
%-------------------------------------------------------------------------------
if nargin<3, dirs=''; end
if nargin<5, exlon=[]; exlat=[]; end
if nargin<6, minland=10; end
if nargin<7, minlake=10; end
if nargin<8, delriver=0; end
if nargin<9, mindp=0; end

src=fullfile(dirs,[coast,'.b']);

f=fopen(src,'r','ieee-be');
if f<0, disp(['input coast line file open error : ',src]), return, end

p=[]; n=0;
while 1
    % read polygon header
    % h(1)   : polygon id number
    % h(2)   : number of points
    % h(3)   : level (1:land,2:lake,3:island_in_lake,...)
    % h(4:7) : min/max extent west/east/south/north (udeg)
    % h(8)   : area of polygon (1/10km^2)
    % h(9)   : greenwich is crossed/source
    [h,cnt]=fread(f,9,'int32'); if cnt<=0, break, end
    
    % read polygon points(lon/lat(udeg))
    pp=fread(f,h(2)*2,'int32');
    plon=pp(1:2:end)*1E-6; i=find(plon>180); plon(i)=plon(i)-360;
    plat=pp(2:2:end)*1E-6;
    pinf=[h(1:3);h(4:7)*1E-6;h(8)/10;floor(h(9)/65536);rem(h(9),65536)];
    
    add=(pinf(3)==1&pinf(8)>=minland)|(pinf(3)>1&pinf(8)>=minlake);
    if add&delriver&pinf(3)>1
        add=(pinf(5)-pinf(4))*(pinf(7)-pinf(6))*1000<pinf(8);
    end
    if add&~isempty(exlat)
        add=exlon(1)<=pinf(5)&pinf(4)<=exlon(2)&exlat(1)<=pinf(7)&pinf(6)<=exlat(2);
    end
    % add polygon
    s=' ';
    if add
        if plon(1)~=plon(end)|plat(1)~=plat(end)
            plon=[plon;plon(1)];
            plat=[plat;plat(1)];
        end
        if mindp>0, [plon,plat,pinf(2)]=reducepoly(plon,plat,mindp,pinf(2)); end
        if length(plon)>2
            n=n+1;
            p(n).lon=plon;
            p(n).lat=plat;
            p(n).inf=pinf;
            s='*';
        end
    end
    msg=sprintf('%s%4d: np=%4d level=%d lon=%5.1f,%5.1f lat=%5.1f,%5.1f area=%8.0f inf=%d,%d',s,pinf);
    disp(msg);
end
fclose(f);

save(file,'p','src','exlon','exlat');

disp(sprintf('generated : %s, no of segment : %d',file,length(p)))

% reduce polygon points ------------------------------------------------------------
function [lon,lat,np]=reducepoly(lon,lat,mindp,np)
dp2=(lon(2:end)-lon(1:end-1)).^2+(lat(2:end)-lat(1:end-1)).^2;
i=dp2<mindp^2;
if length(find(i))>1
    i=find(i(1:end-1)&i(2:end))+1;
    lat(i)=[]; lon(i)=[];
    disp(sprintf('polygon reduction : np=%d->%d',np,length(lat)))
end
