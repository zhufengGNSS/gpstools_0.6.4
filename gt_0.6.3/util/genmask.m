function genmask
% generate japan map mask

db='gshhs_i.mat';
xl=[128,146];
yl=[28,48];

load(db); info=[p.inf]';
n=0; m=0;
for i=find(xl(1)<=info(:,4)&info(:,5)<=xl(2)&yl(1)<=info(:,6)&info(:,7)<=yl(2))'
    if mod(info(i,3),2)==1
        n=n+1; lonc{n}=flipud(p(i).lon); latc{n}=flipud(p(i).lat);
    else
        m=m+1; lonl{m}=flipud(p(i).lon); latl{m}=flipud(p(i).lat);
    end
end
[xs,ys]=mapmask(lonc,latc,xl,yl);
figure, hold on, grid on
patch(xs,ys,'w','edgecolor','none');
for n=1:length(lonl), patch(lonl{n},latl{n},'w','edgecolor','none'); end
xlim(xl); ylim(yl); set(gca,'color','none');

% generate map mask ------------------------------------------------------------
function [xs,ys]=mapmask(x,y,xl,yl)

x={x{:},xl'}; y={y{:},yl([1,1])'}; % add bottom line

for n=1:length(y)-1
    [miny,i]=min(y{n}); % bottom point of map
    mm=[]; jj=[]; pp=[];
    for m=n+1:length(y)
        [j,p]=crossp(x{m},y{m},x{n}(i),y{n}(i));
        if ~isempty(j), mm=[mm;m]; jj=[jj;j]; pp=[pp;p]; end
    end
    [maxy,k]=max(pp(:,2)); m=mm(k); j=jj(k); p=pp(k,:);
    [x{m},y{m}]=catmap(x{n},y{n},x{m},y{m},i,j,p);
end
xs=[x{end};xl([2,1,1])'];
ys=[y{end};yl([2,2,1])'];

% crossing point with lower vertical half-line ---------------------------------
function [i,p]=crossp(x,y,xx,yy)
i=[]; p=[]; if xx<min(x)|max(x)<xx, return, end
for j=find(x(1:end-1)<=xx&xx<x(2:end))' % x(i)<=xx<x(i+1)
    i=[i;j];
    a=(xx-x(j))/(x(j+1)-x(j));
    p=[p;x(j)*(1-a)+x(j+1)*a,y(j)*(1-a)+y(j+1)*a];
end
j=find(p(:,2)<yy); [ymax,k]=max(p(j,2));
i=i(j(k)); p=p(j(k),:);

% connect maps -----------------------------------------------------------------
function [x,y]=catmap(x1,y1,x2,y2,i,j,p)
x=[x2(1:j);p(1);x1(i:end);x1(1:i);p(1);x2(j+1:end)];
y=[y2(1:j);p(2);y1(i:end);y1(1:i);p(2);y2(j+1:end)];
