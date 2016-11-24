function [zpds,sigs]=readtrop(td,time,rcvs,tropdir,tropsrc,tunit,opts)
%-------------------------------------------------------------------------------
% [system] : GpsTools
% [module] : read tropospheric parameters
% [func]   : read tropospheric parameters
% [argin]  : td      = date (mjd-gpst)
%            time    = time vector (sec)
%            rcvs    = station list
%           (tropdir) = data directory (default:current)
%           (tropsrc) = data source    (default:'igssnx')
%                    'igssnx' : IGS sinex(trop)
%                    'igsmon' : IGS monthly sinex(trop)
%                    'zpdf'   : estimated zpd (forward)
%                    'zpdb'   : estimated zpd (backward)
%                    'zpdfb'  : estimated zpd (smoothed)
%                    'mso'    : gpv derived zpd
%                    'saast'  : saastamoinen+standard atmosphere model
%           (tunit)  = processint unit time (H) for 'zpdf','zpdb','zpdfb'
%           (opts)   = options
%                   'interp' : interporate values
% [argout] : zpds = zenith total delay (m)
%                  zpds(n,:,m)=time(n),rcvs{n}zpd
%            dsxs = horizontal correction [dstx,dsty,dstz] (m)
% [note]   :
% [version]: $Revision: 16 $ $Date: 2008-12-12 15:49:30 +0900 (é‡‘, 12 12 2008) $
%            Copyright(c) 2004-2008 by T.Takasu, all rights reserved
% [history]: 04/11/19  0.1  new
%            04/11/22  0.2  add option 'intp'
%            06/02/07  0.3  add 'ztd*','mso' as tropsrc, add argin tunit
%            08/11/24  0.4  fix bug on error-stop if format error
%-------------------------------------------------------------------------------
if nargin<4, tropdir=''; end
if nargin<5, tropsrc=''; end
if nargin<6, tunit=24; end
if nargin<7, opts=''; end
if isempty(tropsrc), tropsrc='igssnx'; end
if ischar(rcvs), rcvs={rcvs}; end
rcvs=upper(rcvs);
zpds=repmat(nan,[length(time),1,length(rcvs)]); sigs=zpds;

switch tropsrc
case {'igssnx','igsmon'}
    for n=1:length(rcvs)
        for m=1:length(time)
            [zpds(m,:,n),sigs(m,:,n)]=...
                ReadTropData(td,time(m),rcvs{n},tropdir,tropsrc,opts);
        end
    end
case {'zpdf','zpdb','zpdfb'}
    zpds=repmat(nan,[length(time),1,length(rcvs)]); sigs=zpds;
    for n=1:length(rcvs)
        [t,xs,vs,prm]=readest(td,time,'zpd',rcvs{n},tropdir,tropsrc(4:end),tunit);
        [tt,i,j]=intersect(time,t);
        if ~isempty(i)
            zpds(i,1,n)=xs(j,1); sigs(i,1,n)=sqrt(vs(j,1));
            if size(xs,2)==3
                zpds(i,2:3,n)=xs(j,2:3); sigs(i,2:3,n)=sqrt(vs(j,2:3));
            elseif size(xs,2)>=4
                zpds(i,2:3,n)=xs(j,3:4); sigs(i,2:3,n)=sqrt(vs(j,3:4));
            end
            if strcmp(opts,'interp')
                i=find(all(~isnan(zpds(:,1,n)),2));
                if ~isempty(i)
                    zpds(:,:,n)=interp1(time(i),zpds(i,:,n),time,'linear','extrap');
                    sigs(:,:,n)=interp1(time(i),sigs(i,:,n),time,'linear','extrap');
                end
            end
        end
    end
case {'rsm','msm','gsm','rso','mso','gso'}
    zpds=GpvZpd(td,time,rcvs,tropdir,tropsrc,opts);
case 'saast'
    for n=1:length(rcvs)
        gpos=eceftogeod(readpos(0,0,rcvs{n},'','approx')');
        zpds(:,1,n)=trop_saast(0,[0,pi/2],gpos);
    end
otherwise
    disp(['warning : trop-delay data source error : ',tropsrc])
end

% read tropospheric data -------------------------------------------------------
function [zpd,sig]=ReadTropData(td,ts,rcv,tropdir,tropsrc,opts)
persistent file tm zpds sigs, if isempty(file), file=''; end
zpd=nan; sig=nan;

t=td+ts/86400; dt=mjdtocal(t);
gpsd=t-44244; gpsw=floor(gpsd/7);
if strcmp(tropsrc,'igssnx')
    f=gfilepath(tropdir,sprintf('%s%04d.zpd',rcv,gpsw),dt,rcv);
else
    f=gfilepath(tropdir,sprintf('%s_%04d_%02d.zpd',rcv,dt(1:2)),dt,rcv);
end
if ~strcmp(f,file)
    file=f;
    [tm,zpds,sigs]=ReadSinexTrop(file,rcv);
    if isempty(tm), return, end
    tm=[0;tm;99999];
    zpds=zpds([1,1:end,end],:);
    sigs=sigs([1,1:end,end],:);
end
if ~isempty(tm)
    if strcmp(opts,'interp')
        zpd=interp1(tm,zpds,t,'linear','extrap');
        sig=interp1(tm,sigs,t,'linear','extrap');
    else
        i=find(abs(tm-t)<1/86400);
        if ~isempty(i), zpd=zpds(i); sig=sigs(i,:); end
    end
end

% read sinex trop file ---------------------------------------------------------
function [time,zpds,sigs]=ReadSinexTrop(file,rcv)
fd=fopen(file,'rt');
if fd<0
    time=[]; zpds=[]; sigs=[];
    disp(['warning : sinex trop file open error : ',file]), return
end
time=zeros(1000,1); zpds=zeros(1000,1); sigs=zeros(1000,1); nt=0; start=0;
while 1
    str=fgetl(fd); if ~isstr(str), break, end
    if findstr(str,'+TROP/SOLUTION'), start=1;
    elseif findstr(str,'-TROP/SOLUTION'), break; end
    if start&str(1)==' '
        if strcmp(SubStr(str,2,4),rcv)
            nt=nt+1;
            time(nt,1)=TimeToMjd(SubStr(str,7,12));
            zpds(nt,1)=StrToNum(str,19,7)*1E-3;
            sigs(nt,1)=StrToNum(str,27,4)*1E-3;
        end
    end
end
fclose(fd);
time=time(1:nt,:); zpds=zpds(1:nt,:); sigs=sigs(1:nt,:);

% gpv derived zpd --------------------------------------------------------------
function zpds=GpvZpd(td,time,rcvs,dirs,src,opts)
zpds=repmat(nan,[length(time),1,length(rcvs)]);
if ~strcmp(src,'mso'), return, end
poss=readpos(td,0,rcvs,'','approx');
for n=1:length(rcvs)
    gpos(n,:)=eceftogeod(poss(1,:,n)');
    gpos(n,3)=geodh(gpos(n,:));
end
ts=(floor(time(1)/10800):ceil(time(end)/10800))*10800;
zs=repmat(nan,[length(ts),1,length(rcvs)]);
for n=1:length(ts)
    for m=1:length(rcvs)
        [pres,temp,humi]=readmet(td,time(n),gpos(m,:),dirs,src);
        if ~isempty(pres), zs(n,1,m)=mettoztd(gpos(m,:),pres,temp,humi); end
    end
end
[t,i,j]=intersect(time,ts);
for n=1:length(rcvs)
    if strcmp(opts,'interp')
        zpds(:,1,n)=interp1(ts,zs(:,1,n),time,'linear','extrap');
    else
        zpds(i,1,n)=zs(j,1,n);
    end
end

% met parameters to ztd --------------------------------------------------------
function ztd=mettoztd(gpos,pres,temp,humi)
Rv=461;      % J/kg/K
k3=3.754E5;  % k3 (K^2/hPa)
kd=23.7146;  % k2-k1*mv/md(k1=77.6,k2=71.98,mv=18.0152,md=28.9644)

% pres(1)=pmsl, temp(1)=tempsurf, humi(1)=humisurf
p=pres(1)*(1-0.0065*gpos(3)/(temp(1)-0.0065*gpos(3)+273.15))^5.257;
zhd=0.0022768*p/(1-0.00266*cos(2*gpos(1)*pi/180)-2.8E-7*gpos(3));
i=flipud(find(p>pres(2:end))+1);
pwv=gpvpwv([p;pres(i)],temp([1;i]),humi([1;i]));
tm=70.2+0.72*(temp(1)+273.15);
zwd=pwv*Rv*(kd+k3/tm)/1E5;
ztd=zhd+zwd;

% read met parameters ----------------------------------------------------------
function [pres,temp,humi]=readmet(td,time,gpos,dirs,gpvsrc)
persistent tdd t1 t2 pm1 pm2 tm1 tm2 hu1 hu2 gpp1 gpp2 gpt1 gpt2 gph1 gph2 p1 p2
if isempty(tdd)|isempty(t1)|td~=tdd|time<t1|t2<=time
    t1=floor(time/3600)*3600;
    if t1==t2
        pm1=pm2; gpp1=gpp2; tm1=tm2; gpt1=gpt2; hu1=hu2; gpth=gph2;
    else
        ts1=floor(t1/3600/6)*6; ft1=rem(t1/3600,6);
        [pm1,gpp1]=readgpv(td,ts1*3600,'pmsl',dirs,gpvsrc,0,ft1);
        [tm1,gpt1]=readgpv(td,ts1*3600,'temp',dirs,gpvsrc,nan,ft1);
        [hu1,gph1]=readgpv(td,ts1*3600,'humi',dirs,gpvsrc,nan,ft1);
    end
    t2=t1+3600*3; ts2=floor(t2/3600/6)*6; ft2=rem(t2/3600,6);
    [pm2,gpp2]=readgpv(td,ts2*3600,'pmsl',dirs,gpvsrc,0,ft2);
    [tm2,gpt2,p1]=readgpv(td,ts2*3600,'temp',dirs,gpvsrc,nan,ft2);
    [hu2,gph2,p2]=readgpv(td,ts2*3600,'humi',dirs,gpvsrc,nan,ft2);
    tdd=td;
end
% interpolate by position
pn1=nan; tn1=nan; hn1=nan; pn2=nan; tn2=nan; hn2=nan;
if ~isempty(gpp1)&~isempty(gpt1)&~isempty(gph1)
    [x1,y1]=gmt('lltogrid',gpos(2),gpos(1),gpp1);
    [x2,y2]=gmt('lltogrid',gpos(2),gpos(1),gpt1);
    [x3,y3]=gmt('lltogrid',gpos(2),gpos(1),gph1);
    pn1=interp2(1:gpp1.nx,1:gpp1.ny,double(pm1),x1,y1);
    for n=1:length(p1)
        tn1(n,1)=interp2(1:gpt1.nx,1:gpt1.ny,double(tm1(:,:,n)),x2,y2);
    end
    for n=1:length(p2)
        hn1(n,1)=interp2(1:gph1.nx,1:gph1.ny,double(hu1(:,:,n)),x2,y2);
    end
end
if ~isempty(gpp1)&~isempty(gpt1)&~isempty(gph2)
    [x1,y1]=gmt('lltogrid',gpos(2),gpos(1),gpp2);
    [x2,y2]=gmt('lltogrid',gpos(2),gpos(1),gpt2);
    [x3,y3]=gmt('lltogrid',gpos(2),gpos(1),gph2);
    pn2=interp2(1:gpp2.nx,1:gpp2.ny,double(pm2),x1,y1);
    for n=1:length(p1)
        tn2(n,1)=interp2(1:gpt2.nx,1:gpt2.ny,double(tm2(:,:,n)),x2,y2);
    end
    for n=1:length(p2)
        hn2(n,1)=interp2(1:gph2.nx,1:gph2.ny,double(hu2(:,:,n)),x2,y2);
    end
end
% interpolate by time
a=(time-t1)/(t2-t1);
if ~isnan(pn1)&~isnan(pn2)
    pmsl=pn1*(1-a)+pn2*a;
    temp=tn1*(1-a)+tn2*a;
    humi=hn1*(1-a)+hn2*a;
elseif a==0&~isnan(pn1)
    pmsl=pn1;
    temp=tn1;
    humi=hn1;
else
    pmsl=nan;
    temp=repmat(nan,length(p1),1);
    humi=repmat(nan,length(p2),1);
end
[pres,i,j]=intersect(p1,p2);
pres(1)=pmsl; temp=temp(i); humi=humi(j);

% compute pwv ------------------------------------------------------------------
function pwv=gpvpwv(pres,temp,humi)
pwv=0; g=9.81;
for n=1:length(pres)-1
    r1=mrat(pres(n),temp(n),humi(n));
    r2=mrat(pres(n+1),temp(n+1),humi(n+1));
    pwv=pwv+1/g*(r1+r2)/2*(pres(n)-pres(n+1))*0.1; % pwv=1/gçrdp
end
if pwv==0, pwv=nan; end

% mixing ratio -----------------------------------------------------------------
function r=mrat(pres,temp,humi)
e=6.11*10^(7.5*temp/(237.3+temp))*humi*0.01; % wv pressure
r=0.622*e/(pres-e);

% substring --------------------------------------------------------------------
function sstr=SubStr(str,pos,ns), sstr=str(pos:min(pos+ns-1,length(str)));

% string->number ---------------------------------------------------------------
function num=StrToNum(str,pos,ns)
num=sscanf(str(pos:min(pos+ns-1,length(str))),'%f',1);
if isempty(num), num=nan; end

% sinex time->mjd---------------------------------------------------------------
function mjd=TimeToMjd(str)
mjd=0; t=sscanf(str,'%f:%f:%f');
if length(t)>=3
    if t(1)<=50, t(1)=t(1)+2000; else t(1)=t(1)+1900; end
    mjd=caltomjd([t(1),1,1])+t(2)-1+t(3)/86400;
end
