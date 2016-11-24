function pwvsonde(year,month,dirs)
%-------------------------------------------------------------------------------
% [system] : GpsTools
% [module] : generate sonde pwv
% [func]   : generate sonde pwv
% [argin]  : year  = year
%            month = months
%            dirs  = sonde data directory
% [argout] : none
% [note]   : <output file>
%            pwvsonde_{rcv}.mat
%                time = time (mjd-jst)
%                data = pwv data
%                    data(n,1)=time(n) pressure (hPa)
%                    data(n,2)=time(n) temperature (degC)
%                    data(n,3)=time(n) PWV (m)
% [version]: $Revision: 12 $ $Date: 2008-11-25 10:02:15 +0900 (ÁÅ´, 25 11 2008) $
%            Copyright(c) 2004-2006 by T.Takasu, all rights reserved
% [history]: 05/06/15   0.1  new
%-------------------------------------------------------------------------------
if nargin<1, year=2004; end
if nargin<2, month=1; end
if nargin<3, dirs=''; end

td=caltomjd([year,month(1),1]);
te=caltomjd([year,month(end)+1,1]);
time=(td:0.5:te-0.5)';

zd=[];
for n=1:length(month)
    file=fullfile(dirs,sprintf('ks%04d%02d.spl',year,month(n)));
    zd=[zd;readsonde(file)];
end
rcvs=unique(zd(:,2));
for n=1:length(rcvs)
    d=zd(zd(:,2)==rcvs(n),:);
    data=repmat(nan,length(time),3);
    for m=1:length(time)
        data(m,:)=genpwv(d(d(:,1)==time(m),3:end));
    end
    save(sprintf('pwvsonde_%03d.mat',rcvs(n)),'time','data');
end

% compute pwv ------------------------------------------------------------------
function d=genpwv(data)
data=data(~isnan(data(:,end)),:);
pwv=0; g=9.81;
for n=1:size(data,1)-1
    r1=mrat(data(n,2),data(n,3),data(n,4));
    r2=mrat(data(n+1,2),data(n+1,3),data(n+1,4));
    pwv=pwv+1/g*(r1+r2)/2*(data(n,2)-data(n+1,2))*0.1; % pwv=1/gÅÁrdp
end
if pwv==0, pwv=nan; end
if isempty(data), d=[nan,nan,pwv]; else d=[data(1,2:3),pwv]; end

% mixing ratio -----------------------------------------------------------------
function r=mrat(pres,temp,humi)
e=6.11*10^(7.5*temp/(237.3+temp))*humi*0.01; % wv pressure
r=0.622*e/(pres-e);

% read sonde file  -------------------------------------------------------------
function data=readsonde(file)
f=fopen(file,'r');
if f<0, disp(['file open error : ',file]); return, end
data=[];
while 1
    d=readrecord(f); if isempty(d), break, end
    data=[data;d];
end
fclose(f);

% read record -------------------------------------------------------------------
function data=readrecord(f)
ps=[1000,925,900,850,800,700,600,500,400,350,300,250,200,175,150,125,100,...
    70,50,40,30,20,15,10,5];
data=[];
while 1
    [a,c]=fread(f,256,'short'); if c<256, break, end
    rcv=double(a(1));
    day=double(a(2:4)); mon=floor(day(2)/100);
    time=caltomjd([day(1),mon,day(2)-mon*100,day(3)-9,0,0]);
    d(1,:)=[time,rcv,0,double(a(7:8))',double(a(9:11))'];
    for n=2:26
        d(n,:)=[time,rcv,double(a(n*7)),ps(n-1)*10,double(a(n*7+1:n*7+4))'];
    end
    d(d<=-32766)=nan;
    i=find(~isnan(d(:,2)));
    
    % [time,rcv,h,pres,temp,humi]
    data=[data;d(i,1:3),d(i,4:5)/10,d(i,6)];
end
