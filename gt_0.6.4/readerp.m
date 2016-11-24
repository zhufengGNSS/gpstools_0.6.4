function erp_value=readerp(tutc,erpdir,erpsrc,tunit)
%-------------------------------------------------------------------------------
% [system] : GpsTools
% [module] : read earth rotation parameters
% [func]   : read earth rotation parameters
% [argin]  : tutc = time(mjd-utc)
%            (erpdir)  = erp directory (default:current)
%            (erpsrc)  = erp source (default:'bulb')
%                        'bula' = IERS Bul.A
%                        'bulb' = IERS Bul.B
%                        'igs'  = IGS Final
%                        'igr'  = IGS Rapid
%                        'c04'  = IERS C04 series
%                        'erpf' = estimated (forward)
%                        'erpb' = estimated (backward)
%                        'erpfb'= estimated (smooothd)
%                        'cod','jpl',... = CODE,JPL,...
%            (tunit)   = processint unit time (hr)
% [argout] : erp_value = earth rotation parameter values
%                 erp_value(n,1) = time(n) pole offset Xp (rad)
%                 erp_value(n,2) = time(n) pole offset Yp (rad)
%                 erp_value(n,3) = time(n) UT1-UTC (sec)
%                 erp_value(n,4) = time(n) nutation offset dpsi (rad)
%                 erp_value(n,5) = time(n) nutation offset deps (rad)
% [note]   : not include daily,sub-daily variations
% [version]: $Revision: 16 $ $Date: 2008-12-12 15:49:30 +0900 (金, 12 12 2008) $
%            Copyright(c) 2004-2006 by T.Takasu, all rights reserved
% [history]: 04/11/02   0.1  new
%            05/01/18   0.2  read previous month bul.B if no bul.B file
%            05/05/31   0.3  add erpsrc 'c04' IERS C04 series
%                            add erpsrc erpf,erpb,erpfb
%            05/09/26   0.4  support cod,jpl,... as erp source
%            06/02/10   0.5  add extrapolation of erp
%            06/06/24   0.6  add warning messages
%            08/11/26   0.7  fix bug on error-stop if format error
%-------------------------------------------------------------------------------
persistent sec2 ts te
if nargin<2, erpdir=''; end
if nargin<3, erpsrc=''; end
if nargin<4, tunit=24; end
if isempty(erpsrc), erpsrc='bulb'; end
utc_tai=prm_utc_tai(tutc);
erp_value=zeros(1,5);

switch erpsrc
case 'bula'
    gt_log('erp type not supported  : %s',erpsrc);

case 'bulb'
    if isempty(sec2)|tutc<ts|te<=tutc, [sec2,ts,te]=ReadBulb(tutc,erpdir); end
    if ~isempty(sec2), erp_value=BulbToErp(tutc,sec2); end

case 'c04'
    erp_value=ReadIersC04(tutc,erpdir,utc_tai);

case {'erpf','erpb','erpfb'}
    erp_value=ReadEstErp(tutc,erpdir,erpsrc,utc_tai,tunit);

otherwise
    erp_value=ReadIgsEop(tutc,erpsrc,erpdir,utc_tai);

end

% read iers bull.b -------------------------------------------------------------
function [sec2,ts,te]=ReadBulb(tutc,erpdir)
persistent file
sec2=[]; ts=[]; te=[];
dt=mjdtocal(tutc-4);
mon=(dt(1)-1988)*12+dt(2);
f=gfilepath(erpdir,sprintf('bulletinb.%03d',mon),dt);
fd=fopen(f,'rt');
if fd<0
    dt(2)=dt(2)-1; if dt(2)<=0, dt(1:2)=[dt(1)-1,12]; end
    f=gfilepath(erpdir,sprintf('bulletinb.%03d',mon-1),dt);
    fd=fopen(f,'rt');
    if fd<0
        if ~strcmp(f,file), gt_log('no earth rotation prms  : %s',f); end
        file=f; return
    end
end
while 1
    buff=fgetl(fd);
    if ~isstr(buff)|findstr(buff,'2 - SMOOTHED VALUES'), break, end
end
n=0;
while 1
    buff=fgetl(fd);
    if ~isstr(buff)|findstr(buff,'3 - NORMAL VALUES'), sec2=sec2(1:n,:); break, end
    p=str2num(buff(10:end));
    if length(p)==8, n=n+1; sec2(n,:)=p; end
end
if n>=1
    ts=sec2(1,1);
    te=sec2(n,1);
    sec2=sec2(1:n,:);
end
fclose(fd);

% iers bul.b->erp --------------------------------------------------------------
function erp_value=BulbToErp(tutc,sec2)
sec2rad=4.8481368110953598E-6;
if tutc<sec2(1,1)|sec2(end,1)<tutc, erp_value=[]; return, end
erps=sec2(:,[1:4,7:end]);
erp=InterpL(erps(:,1),erps(:,2:end),tutc);
erp_value=zeros(1,5);
erp_value(1)=erp(1)*sec2rad; % xp(rad)
erp_value(2)=erp(2)*sec2rad; % yp(rad)
erp_value(3)=erp(3);         % UT1-UTC(sec)
erp_value(4)=erp(4)*1E-3*sec2rad; % dpsi(rad)
erp_value(5)=erp(5)*1E-3*sec2rad; % deps(rad)

% Lagrange interpolation (4points/3degrees) ------------------------------------
function yi=InterpL(x,y,xi)
n=floor((xi-x(1))/(x(2)-x(1)))+1;
n0=max(n-1,1);
n1=min(n+2,size(x,1));
yi=zeros(size(y,2),1);
for n=n0:n1
    nn=n0:n1; nn(n-n0+1)=[];
    yi=yi+prod((xi-x(nn))./(x(n)-x(nn))).*y(n,:)';
end 

% read igs erp -----------------------------------------------------------------
function erp_value=ReadIgsEop(tutc,erpsrc,erpdir,utc_tai)
persistent file1 file2 mjd1 erp1 derp1 mjd2 erp2 derp2
tgps=tutc-(19+utc_tai)/86400; dt1=mjdtocal(tgps); dt2=mjdtocal(tgps+7);
if isempty(file1), file1=''; end, if isempty(file2), file2=''; end
erp_value=zeros(1,5);
if strcmp(erpsrc,'igr')
    [gpsw1,gpsd1]=MjdToGpsD(tgps-0.5); dt1=mjdtocal(tgps-0.5);
    [gpsw2,gpsd2]=MjdToGpsD(tgps+0.5); dt2=mjdtocal(tgps+0.5);
    f1=sprintf('igr%04d%d.erp',gpsw1,gpsd1);
    f2=sprintf('igr%04d%d.erp',gpsw2,gpsd2);
    if ~strcmp(f1,file1), [mjd1,erp1,derp1]=ReadIgsEopFile(erpdir,f1,dt1); file1=f1; end
    if ~strcmp(f2,file2), [mjd2,erp2,derp2]=ReadIgsEopFile(erpdir,f2,dt2); file2=f2; end
    if ~isempty(mjd1)&~isempty(mjd2)
        erp_value=erp1*(mjd2-tgps)+erp2*(tgps-mjd1);
    elseif ~isempty(mjd1)
        erp_value=erp1+derp1*(tgps-mjd1);
    elseif ~isempty(mjd2)
        erp_value=erp2+derp2*(tgps-mjd2);
    end
else
    [gpsw,gpsd]=MjdToGpsD(tgps-0.5);
    f1=sprintf('%s%04d7.erp',erpsrc,gpsw);
    f2=sprintf('%s%04d7.erp',erpsrc,gpsw+1);
    if ~strcmp(f1,file1), [mjd1,erp1,derp1]=ReadIgsEopFile(erpdir,f1,dt1); file1=f1; end
    if ~strcmp(f2,file2), [mjd2,erp2,derp2]=ReadIgsEopFile(erpdir,f2,dt2); file2=f2; end
    if isempty(mjd1)&isempty(mjd2), return, end
    [m,i]=intersect(mjd1,mjd2); mjd1(i)=[]; erp1(i,:)=[];
    erp_value=interp1([mjd1;mjd2],[erp1;erp2],tgps,'linear','extrap');
end
erp_value=[erp_value(1:2)*1E-6*pi/180/3600,erp_value(3)*1E-7,0,0];

% read iers c04 ----------------------------------------------------------------
function erp_value=ReadIersC04(tutc,erpdir,utc_tai)
persistent file mjd erp, if isempty(file), file=''; mjd=[]; erp=[]; end
erp_value=zeros(1,5);
dt1=mjdtocal(tutc);
dt2=mjdtocal(tutc+0.99999999);
f1=gfilepath(erpdir,sprintf('eopc04.%02d',mod(dt1(1),100)),dt1);
f2=gfilepath(erpdir,sprintf('eopc04.%02d',mod(dt2(1),100)),dt2);
if ~strcmp([f1,f2],file)
    [mjd,erp]=ReadIersC04File(f1); if isempty(mjd), return, end
    if ~strcmp(f1,f2)
        [mjd2,erp2]=ReadIersC04File(f2); mjd=[mjd;mjd2]; erp=[erp;erp2];
    end
    file=[f1,f2];
end
if tutc<mjd(1)|mjd(end)<tutc
    gt_log('erp data out of range   : %s',f1); return
end
erp_value=interp1(mjd,erp,tutc,'linear','extrap');
erp_value=[erp_value(1:2)*pi/180/3600,erp_value(3),erp_value(4:5)*pi/180/3600];

% read igs erp file ------------------------------------------------------------
function [mjd,erp,derp]=ReadIgsEopFile(erpdir,file,dt)
mjd=[]; erp=[]; derp=zeros(1,3);
file=gfilepath(erpdir,file,dt);
fd=fopen(file,'rt');
if fd<0, gt_log('no earth rotation prms  : %s',file); return, end
n=-1;
while 1
    buff=fgetl(fd); if ~isstr(buff), break, end
    if findstr(buff,'  MJD'), n=0;
    elseif n==0, n=n+1;
    elseif n>0
        value=sscanf(buff,'%f');
        if length(value)>=14
            mjd(n,1)=value(1); erp(n,:)=value(2:4)'; derp(n,:)=[value(13:14)',0,];
            n=n+1;
        end
    end
end
fclose(fd);
if isempty(mjd), gt_log('erp read error          : %s',file); end

% read iers c04 file -----------------------------------------------------------
function [mjd,erp]=ReadIersC04File(file)
mjd=[]; erp=[];
fd=fopen(file,'rt');
if fd<0, gt_log('no earth rotation prms  : %s',file); return, end
n=-1;
while 1
    buff=fgetl(fd); if ~isstr(buff), break, end
    if findstr(buff,'YEAR ==>'), n=0;
    elseif n==0, n=n+1;
    elseif n>0&~isempty(deblank(buff))
        mjd(n,1)=StrToNum(buff,12,5);
        erp(n,1)=StrToNum(buff,17,9);
        erp(n,2)=StrToNum(buff,26,9);
        erp(n,3)=StrToNum(buff,35,10);
        erp(n,4)=StrToNum(buff,59,9);
        erp(n,5)=StrToNum(buff,68,9);
        n=n+1;
    end
end
fclose(fd);
if isempty(mjd), gt_log('erp read error          : %s',file); end

% read estimated erp -----------------------------------------------------------
function erp=ReadEstErp(tutc,erpdir,erpsrc,utc_tai,tunit)
persistent ts xs vs dirs src
t=tutc-(19+utc_tai)/86400;
if isempty(ts)|t<ts(1)|ts(end)<t|~strcmp(dirs,erpdir)|~strcmp(src,erpsrc)
   td=floor(t);
   [ts,xs,vs]=readest(td,mod(t,1)*86400,'erp','',erpdir,erpsrc(4:end),tunit);
   ts=td+ts/86400; dirs=erpdir; src=erpsrc;
end
if ~isempty(ts)
    erp=interp1(ts,xs,t,'linear','extrap');
else
    erp=[nan,nan,nan];
end
erp=[erp(1:2)*pi/180/3600,erp(3),0,0];

% MJD->GPS Week/day/hour -------------------------------------------------------
function [gpsw,gpsd,hour]=MjdToGpsD(mjd)
gpsd=mjd-44244; gpsw=floor(gpsd/7); gpsd=gpsd-gpsw*7;
hour=mod(gpsd,1)*24; gpsd=floor(gpsd);

% string->number ---------------------------------------------------------------
function num=StrToNum(str,pos,ns)
num=sscanf(str(pos:min(pos+ns-1,length(str))),'%f',1);
if isempty(num), num=0; end
