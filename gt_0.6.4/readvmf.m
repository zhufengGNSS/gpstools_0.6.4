function [ah,aw,zh,zw]=readvmf(time,gpos,dirs)
%-------------------------------------------------------------------------------
% [system] : GpsTools
% [module] : read VMF1 coefficients
% [func]   : read VMF1 coefficients
% [argin]  : time  = date-time (mjd-utc)
%            gpos  = receiver position [lat,lon,height] (deg,m)
%           (dirs) = VMF1 coefficients directory (default: current)
% [argout] : ah    = hydrostatic coefficient a
%            aw    = wet coefficient a
%           (zh)   = zenith hydrostatic delay (m)
%           (zw)   = zenith wet delay (m)
% [note]   :
% [version]: $Revision:$ $Date:$
%            Copyright(c) 2004-2008 by T.Takasu, all rights reserved
% [history]: 08/12/06   0.1  new
%-------------------------------------------------------------------------------
persistent dirs1_ time1_ ahs1 ahs2 aws1 aws2
persistent dirs2_ time2_ zhs1 zhs2 zws1 zws2
if nargin<3, dirs=''; end

if any(isnan(gpos(1:2)))
    ah=nan; aw=nan; zh=nan; zw=nan;
    return;
end
% read vmf coefficient files
t=floor(time/0.25)*0.25; tt=(time-t)/0.25;
if isempty(time1_)|t~=time1_|~strcmp(dirs,dirs1_)
    file1=vmffile(t); file2=vmffile(t+0.25);
    ahs1=readvmfgrid(fullfile(dirs,['ah',file1]));
    aws1=readvmfgrid(fullfile(dirs,['aw',file1]));
    if isempty(ahs1)|isempty(aws1)
        gt_log('no vmf1 grid files      : %s',fullfile(dirs,['ah|aw',file1]));
    end
    ahs2=readvmfgrid(fullfile(dirs,['ah',file2]));
    aws2=readvmfgrid(fullfile(dirs,['aw',file2]));
    if isempty(ahs2)|isempty(aws2)
        gt_log('no vmf1 grid files      : %s',fullfile(dirs,['ah|aw',file2]));
    end
    dirs1_=dirs; time1_=t;
end
if nargout>=3&(isempty(time2_)|t~=time2_)|~strcmp(dirs,dirs2_)
    file1=vmffile(t); file2=vmffile(t+0.25);
    zhs1=readvmfgrid(fullfile(dirs,['zh',file1]));
    zws1=readvmfgrid(fullfile(dirs,['zw',file1]));
    zhs2=readvmfgrid(fullfile(dirs,['zh',file2]));
    zws2=readvmfgrid(fullfile(dirs,['zw',file2]));
    dirs2_=dirs; time2_=t;
end
ah1=interppos(gpos,ahs1);
ah2=interppos(gpos,ahs2);
aw1=interppos(gpos,aws1);
aw2=interppos(gpos,aws2);
ah=interptime(tt,ah1,ah2)/1E8;
aw=interptime(tt,aw1,aw2)/1E8;
if nargout>=3
    zh1=interppos(gpos,zhs1);
    zh2=interppos(gpos,zhs2);
    zw1=interppos(gpos,zws1);
    zw2=interppos(gpos,zws2);
    zh=interptime(tt,zh1,zh2);
    zw=interptime(tt,zw1,zw2);
end

% vmf grid file name -----------------------------------------------------------
function file=vmffile(time)
dv=mjdtocal(time);
day=time-caltomjd([dv(1),1,1])+1;
doy=floor(day); hr=(day-doy)*24;
file=sprintf('%02d%03d.h%02.0f',mod(dv(1),100),doy,hr);

% interpolate coefficients by position -----------------------------------------
function c=interppos(gpos,cs)
if gpos(2)<0, gpos(2)=gpos(2)+360; end
if isempty(cs)|gpos(1)<-90|90<gpos(1)|gpos(2)<0|360<=gpos(2), c=nan; return; end
% [lat,lon]=[90:-2:-90,0:2.5:360]
a=(90-gpos(1))/2; i=min(floor(a),89); a=a-i;
b=gpos(2)/2.5; j=floor(b); b=b-j;
c=cs(i+1,j+1)*(1-a)*(1-b)+cs(i+2,j+1)*a*(1-b)+cs(i+1,j+2)*(1-a)*b+cs(i+2,j+2)*a*b;

% interpolate coefficients by time ---------------------------------------------
function c=interptime(tt,c1,c2)
if ~isnan(c1)&~isnan(c2), c=c1*(1-tt)+c2*tt;
elseif ~isnan(c1), c=c1;
elseif ~isnan(c2), c=c2;
else c=nan; end
