function testeph(td,varargin)
if nargin<1, td=caltomjd([2003,12,1]); end
eph='ephs'; ref='igs'; xyz=0; plotf=0; n=1;
while n<=length(varargin)
    switch varargin{n}
    case 'eph', eph=varargin{n+1}; n=n+2;
    case 'ref', ref=varargin{n+1}; n=n+2;
    case 'xyz', xyz=1; n=n+1;
    case 'plot', plotf=1; n=n+1;
    otherwise, error('argin error'), end
end
for tdd=td, teste(tdd,eph,ref,xyz,plotf), end

function teste(td,eph,ref,xyz,plotf)
gpsd=td-44244; gpsw=floor(gpsd/7); gpsd=gpsd-gpsw*7;
C=299792458;
file=sprintf('%4d%d',gpsw,gpsd);
[epr,timer,ephr,satr]=readsp3(['d:\gps\eph\',ref,file,'.sp3']);
[eps,times,ephs,sats]=readsp3([eph,file,'.sp3']);
if isempty(eps)
    [eps,times,ephs,sats]=readsp3([eph,file,'.eph']);
end
if isempty(eps)
    [eps,times,ephs,sats]=readsp3([eph,file,'.sp3c']);
end
if isempty(eps)
    [eps,times,ephs,sats]=readsp3([eph,file,'_00.sp3']);
end
ep=mjdtocal(td);
day=sprintf('%d/%d/%d',ep(1:3));
merr=repmat(nan,length(sats),4);

if xyz, label={'3d(m)','x(m)','y(m)','z(m)'};
else label={'3d(m)','radial(m)','alongt(m)','crosst(m)'}; end

if ~plotf, msg('%-12s %-7s %8s %8s %8s %8s','   day','  sat',label{:}), end

for n=1:length(sats)
    i=find(strcmp(sats{n},satr));
    if ~isempty(i)
        err=ephs(:,1:3,n)-ephr(:,1:3,i); % ecef
        if size(ephs,2)>=6&~xyz
            for m=1:size(ephs,1)
                [U,P,N,gmst,dx,dy,du]=ecsftoecef(td+(times(m)+19-32)/86400);
                pos=ephs(m,1:3,n)*U;
                vel=ephs(m,4:6,n)*U+ephs(m,1:3,n)*du;
                err(m,:)=err(m,:)*U*ecsftosatf([pos,vel]')';
            end
        end
        err=[sqrt(sum(err(:,1:3).^2,2)),err];
        for m=1:4, merr(n,m)=rmsn(err(:,m)); end
    end
    if ~plotf
        msg('%-12s %-7s %8.4f %8.4f %8.4f %8.4f',day,sats{n},merr(n,:))
    else
        figure
        for m=1:4
            subplotv(m,4,td,[0,24],[-0.25,0.25],label{m})
            if m==1, ylim([0,0.5]), title([sats{n},' ',day]), end
            plot(times/3600,err(:,m))
        end
    end
end
if ~plotf
    msg('%-12s %-7s %8.4f %8.4f %8.4f %8.4f',day,'Average',...
        meann(merr(:,1)),meann(merr(:,2)),meann(merr(:,3)),meann(merr(:,4)))
end

function m=meann(x)
i=find(~isnan(x)); if isempty(i), m=nan; else m=mean(x(i)); end

function m=rmsn(x)
i=find(~isnan(x)); if isempty(i), m=nan; else m=sqrt(mean(x(i).^2)); end

function subplotv(n,m,td,tspan,range,yl)
subplot(m,1,n), hold on, grid on, box on, ylim(range), taxis(td), ylabel(yl)

function msg(varargin), disp(sprintf(varargin{:}))
