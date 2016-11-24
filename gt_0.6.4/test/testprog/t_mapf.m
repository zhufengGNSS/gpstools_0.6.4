function t_mapf(td,time,gpos,azel,dirs)
% test mapping functions
if nargin<1, td=caltomjd([2006,1,1]); end
if nargin<2, time=0:10800:86400*731; end
if nargin<3, gpos=[-36.1,140.1,67]; end
if nargin<4, azel=[0,5]; end
if nargin<5, dirs='g:\vmf'; end

dv=mjdtocal(td);
ty=dv(1)+(td-caltomjd([dv(1),1,1]))/365.25+time/86400/365.25;
azel=azel*pi/180;

mfh=repmat(nan,length(time),3);
mfw=repmat(nan,length(time),3);
for i=1:length(time)
    mjd=td+time(i)/86400; dv=mjdtocal(mjd);
    disp(sprintf('%04d/%02d/%02d %02.0fH',dv(1:3),mod(time(i)/3600,24)));
    %[ah,aw]=readvmf(mjd,gpos,dirs);
    [mfh(i,1),mfw(i,1)]=mapf_nmf(mjd,azel,gpos);
    [mfh(i,2),mfw(i,2)]=mapf_gmf(mjd,azel,gpos);
    %[mfh(i,3),mfw(i,3)]=mapf_vmf1(mjd,azel,gpos,ah,aw);
end
plotmf(ty,mfh,'Hydrostatic Mapping Function',gpos,azel)
plotmf(ty,mfw,'Wet Mapping Function',gpos,azel)

% plot mapping function
function plotmf(ty,mf,ti,gpos,azel)
fname='Times New Roman';
fsize=9;
figure('color','w'); hold on; box on;
plot(ty,mf');
set(gca,'fontname',fname,'fontsize',fsize);
xlabel('Year');
ylabel(sprintf('Mapping Function (El=%.0f\\circ)',azel(2)*180/pi));
xlim(ty([1,end]));
legend({'NMF','GMF','VMF1'});
title(sprintf('%s (Lat: %.1f\\circ Lon: %.1f\\circ Height: %.0fm)',ti,gpos));
moveax
