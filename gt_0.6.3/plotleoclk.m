function plotleoclk(td,time,sat,nsig,dirs)
% plot GRACE orbit error
if nargin<1, td=caltomjd([2004,1,1]); end
if nargin<2, time=0:30:86370; end
if nargin<3, sat='GRACE-A'; end
if nargin<4, nsig=0; end
if nargin<5, dirs=''; end
C=299792458;
dirclk='clk';

[clks,covs]=readclk(td,time,['GRA',sat(end)],dirs,'clkfb',24);
i=find(covs(:,1)<1); time=time(i); clks=clks(i,:); covs=covs(i,:);

t=[]; p=[];
for n=floor(time(1)/86400):floor(time(end)/86400)
    ep=mjdtocal(td+n);
    file=sprintf('CLK1B_%04d-%02d-%02d_%s_00.dat',ep(1:3),sat(end));
    [e,tt,k,pp]=readgrace1b(fullfile(dirclk,file));
    t=[t;tt+(caltomjd(e)-td)*86400];
    p=[p;pp];
end
[tt,i,j]=intersect(time,t);
err=clks(i,1)-p(j,1);
gut('newfigg','','',[600,400]);
margin=[0.08,0.03,0.04,0.012];
ggt('subplotv','',1,2,'taxis',td,'topts','nolabel','xlim',[time(1),time(end)]/3600,...
    'ylim',[-inf,inf],'margin',margin);
plot(tt/3600,clks(i,1)*1E9);
plot(t/3600,p*1E9,'r');

title(['Satellite Clock ',sat,' : ',tstr(td+time(1)/86400),'-',tstr(td+time(end)/86400)],...
      'fontsize',9);
ylabel('Clock Bias (nsec)');

ggt('subplotv','',2,2,'taxis',td,'xlim',[time(1),time(end)]/3600,'ylim',[-0.25,0.25],...
    'margin',margin);
h=plot(tt/3600,err*1E9);
ylabel('Clock Bias Error (nsec)');
bias=mean(err);
rmse=sqrt(mean(err.^2,1));
ggt('mtext','_2',sprintf('REF: Level1B\nMEAN %.4fm RMS %.4fm',bias,rmse),3);
set(gcf,'name',['LEO Satellite Orbit : ',dirs])

function s=tstr(t)
ep=mjdtocal(t); s=sprintf('%04d/%02d/%02d %02d:%02d',ep(1:5));
