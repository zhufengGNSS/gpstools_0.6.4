function metprm=readmet(td,time,rcvs,metdir,metsrc)
%-------------------------------------------------------------------------------
% [system] : GpsTools
% [module] : read meterological parameters
% [func]   : read meterological parameters
% [argin]  : td   = day (mjd-gpst)
%            time = time vector (sec)
%            rcvs = station names
%            (metdir)  = meteorological parameters directory (default:current)
%            (metsrc)  = meteorologibal parameters source (default:'mso')
%                        'mso'= JMA MSM Online
% [argout] : metprm = meterological parameters
%                 metprm(n,1,m) = rcvs{n} time(m) pressure (hPa)
%                 metprm(n,2,m) = rcvs{n} time(m) temperture (C)
%                 metprm(n,3,m) = rcvs{n} time(m) relative humidity (%)
% [note]   :
% [version]: $Revision: 2 $ $Date: 06/07/08 14:21 $
%            Copyright(c) 2004-2006 by T.Takasu, all rights reserved
% [history]: 05/06/15   0.1  new
%-------------------------------------------------------------------------------
if nargin<4, metdir=''; end
if nargin<5, metsrc='mso'; end

switch metsrc
case {'rms','msm','gsm','rso','mso','gso'}
    metprm=ReadGpvData(td,time,rcvs,metdir,metsrc);

otherwise
    disp(['warning : met data source error : ',metsrc])
end

% read gpv --------------------------------------------------------------------
function metprm=ReadGpvData(td,time,rcvs,metdir,metsrc)
metprm=repmat(nan,[length(rcvs),3,length(time)]);
for n=1:length(rcvs)
    gpos(n,:)=eceftogeod(readpos(td,time(1),rcvs{n},'','approx')');
    gpos(n,3)=geodh(gpos(n,:));
end
switch metsrc
case {'mso','gso'}
    if strcmp(metsrc,'mso'), tu=6*3600; ftt=0:5; else tu=12*3600; ftt=0:6:6; end
    m=0; t=(floor(time(1)/tu):floor(time(end)/tu)+1)*tu;
    for n=1:length(t)
        for ft=ftt
            m=m+1; ts(m)=t(n)+ft*3600;
            [pmsl(:,:,m),gprm]=readgpv(td,t(n),'pmsl',metdir,metsrc,0,ft);
            [temp(:,:,m)]=readgpv(td,t(n),'temp',metdir,metsrc,0,ft);
            [humi(:,:,m)]=readgpv(td,t(n),'humi',metdir,metsrc,0,ft);
            if ts(m)>time(end), break, end
        end
    end
    if ~isempty(gprm)
        [x,y,z]=meshgrid(1:gprm.nx,1:gprm.ny,ts);
        pmsl=double(pmsl);
        temp=double(temp);
        humi=double(humi);
        for n=1:length(rcvs)
            [xi,yi]=gmt('lltogrid',gpos(n,2),gpos(n,1),gprm);
            pm=shiftdim(interp3(x,y,z,pmsl,xi,yi,time),2);
            te=shiftdim(interp3(x,y,z,temp,xi,yi,time),2);
            hu=shiftdim(interp3(x,y,z,humi,xi,yi,time),2);
            pr=pm.*(1-0.0065*gpos(n,3)./(te+273.15+0.0065*gpos(n,3))).^5.257;
            metprm(n,:,:)=[pr,te,hu]';
        end
    end
otherwise
    disp(['warning : met data source error : ',metsrc])
end
