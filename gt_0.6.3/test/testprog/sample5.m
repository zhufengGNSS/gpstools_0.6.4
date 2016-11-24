% マルチパスマップ(格子版)
function sample5

dirs.obs='d:\gps\obs\rinex';
dirs.nav='d:\gps\nav';
dirs.pos='d:\gps\pos';

td=CalToMjd([2003,12,1]);
%time=0:30:86400*31;
time=0:30:86400*3;

sats={...
'GPS01','GPS02','GPS03','GPS04','GPS05','GPS06','GPS07','GPS08','GPS09','GPS10',...
'GPS11','GPS13','GPS14','GPS15','GPS16','GPS17','GPS18','GPS20','GPS21','GPS23',...
'GPS24','GPS25','GPS26','GPS27','GPS28','GPS29','GPS30','GPS31'...
};
rcv='GUAM';

es=5;         % 格子間隔(仰角)(度)
as=5;         % 格子間隔(方位角)(度)
narcmax=2000; % アーク最大数
range=[-1,1]; % 表示レンジ

ne=90/es+1; na=360/as;
[z,iz]=ReadObs(td,time,sats,rcv,dirs.obs); if isempty(iz), return, end

AA=zeros(ne*na+narcmax); Ay1=zeros(ne*na+narcmax,1); Ay2=Ay1; narc=0;

for n=1:length(sats)
    i=find(iz(:,2)==n);
    if ~isempty(i)
        [t,zz,arc]=ExtractArc(iz(i,1),z(i,:),120,3,300);
        azels{n}=[];
        for a=arc'
            narc=narc+1;
            [mp1,mp2]=MultiPath(zz(a(1):a(2),:));
            disp(sprintf('%s : arc=%7d-%7dsec : std(mp1)=%6.3fm std(mp2)=%6.3fm',...
                 sats{n},t(a(1)),t(a(2)),std(mp1),std(mp2)))
            azel=SatDir(td,t(a(1):a(2)),sats{n},rcv,dirs);
            sig=1./sin(azel(:,2));
            A=sparse(size(mp1,1),narcmax); A(:,narc)=1;
            A=[DesignMat(azel,na,ne),A];
            A=A./repmat(sig,1,size(A,2));
            AA=AA+A'*A;
            Ay1=Ay1+A'*(mp1./sig);
            Ay2=Ay2+A'*(mp2./sig);
            azels{n}=[azels{n};azel];
        end
    end
end
warning off
x1=repmat(nan,ne*na+narcmax,1); x2=x1;
i=find(any(AA));
x1(i)=AA(i,i)\Ay1(i);
x2(i)=AA(i,i)\Ay2(i);
map1=reshape(x1(1:ne*na),ne,na);
map2=reshape(x2(1:ne*na),ne,na);
warning on

dt1=MjdToCal(td+time(1)/86400); dt2=MjdToCal(td+time(end)/86400);
ti=sprintf('%04d/%02d/%02d %02d:%02d-%02d/%02d %02d:%02d',dt1(1:5),dt2(2:5));

PlotMp(map1,azels,range), title(['CODE MULTI-PATH (C1) (m) : ',ti,' : ',rcv])
PlotMp(map2,azels,range), title(['CODE MULTI-PATH (P2) (m) : ',ti,' : ',rcv])

% 衛星方位仰角 -----------------------------------------------------------------
function azel=SatDir(td,time,sat,rcv,dirs)
nav=ReadNav(td,0:86400:time(end),sat,'',dirs.nav);
posr=ReadPos(td,time(1),rcv,dirs.pos)';
azel=zeros(length(time),2);
for n=1:length(time)
    azel(n,:)=SatAzel(NavToState(td,time(n),nav),posr);
end

% アーク分割 -------------------------------------------------------------------
function [t,z,arc]=ExtractArc(t,z,gapmax,slipmin,arcmin)
C=299792458; f1=1.57542E9; f2=1.22760E9; lam1=C/f1; lam2=C/f2;
i=find(~isnan(z(:,1))); t=t(i); z=z(i,:);
zg=lam1*z(:,1)-lam2*z(:,2);
i=find((t(2:end)-t(1:end-1))>gapmax|abs(zg(2:end)-zg(1:end-1))>slipmin);
arc=[1,i'+1;i',size(t,1)]';
arc=arc(find(t(arc(:,2))-t(arc(:,1))>=arcmin),:);

% マルチパス計算 ---------------------------------------------------------------
function [mp1,mp2]=MultiPath(z)
C=299792458; f1=1.57542E9; f2=1.22760E9; lam1=C/f1; lam2=C/f2;
mp1=repmat(nan,size(z,1),1); mp2=mp1;
for n=1:size(z,1)
    ion=-(lam1*z(n,1)-lam2*z(n,2))/(1-f1^2/f2^2);
    mp1(n)=z(n,3)-lam1*z(n,1)-2*ion;
    mp2(n)=z(n,4)-lam2*z(n,2)-2*f1^2/f2^2*ion;
end
mp1=mp1-mean(mp1(find(~isnan(mp1))));
mp2=mp2-mean(mp2(find(~isnan(mp2))));

% 計画行列 ---------------------------------------------------------------------
function A=DesignMat(azel,na,ne)
A=sparse(size(azel,1),ne*na);
for n=1:size(azel,1)
    A(n,:)=reshape(SphGridFunc(azel(n,:),na,ne),1,ne,na);
end

% 球面グリッド関数(半球) -------------------------------------------------------
function g=SphGridFunc(azel,na,ne)
g=sparse(ne,na); azel=azel*180/pi; es=90/(ne-1); as=360/na;
n=floor(azel(2)/es);
m=floor((azel(1)+180)/as);
a=(azel(2)-es*n)/es;
b=(azel(1)+180-as*m)/as;
gg=[(1-a)*(1-a),a*(1-b);(1-a)*b,a*b];
if m+1<na, g(n+1:n+2,m+1:m+2)=gg; else g(n+1:n+2,[m+1,1])=gg; end

% マップ表示 -------------------------------------------------------------------
function PlotMp(map,azels,range)
es=90/(size(map,1)-1); as=360/size(map,2);
map=[map,map(:,1)];
figure, axis([-180,180,0,90])
for n=1:size(map,1)
    el=(n-1)*es; ry=[el-es/2,el+es/2,el+es/2,el-es/2];
    for m=1:size(map,2)
        az=(m-1)*as-180; rx=[az-as/2,az-as/2,az+as/2,az+as/2];
        if ~isnan(map(n,m)), patch(rx,ry,map(n,m),'edgecolor','none'); end
        line([az,az],[-90,90],'linestyle',':')
    end
    line([-180,180],[el,el],'linestyle',':')
end
set(gca,'xtick',-180:30:180,'ytick',0:15:90)
xlabel('Azimath (deg)'), ylabel('Elevation (deg)')
caxis(range), h=colorbar('horiz');
p=get(h,'position'); set(h,'position',[p(1:3),p(4)*0.3])
hold on
for n=1:length(azels)
    plot(azels{n}(:,1)*180/pi,azels{n}(:,2)*180/pi,'.','markersize',1)
end
