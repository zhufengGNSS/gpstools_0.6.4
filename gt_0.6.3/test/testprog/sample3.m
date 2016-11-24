% TECマップ表示
function sample3

file='d:\gps\ion\igsg3350.03i';

[epoch,time,tec,rms,lats,lons,hgts,rb]=ReadIonex(file);

[td,ts]=CalToMjd(epoch);

for n=1:length(time)
    figure
    contourf(lons,lats,tec(:,:,1,n),0:2.5:100,'b')
    
    caxis([0,100]), h=colorbar('horiz');
    p=get(h,'position'); set(h,'position',[p(1:3),p(4)*0.3])
    
    title(sprintf('IONOSPHERE MAP : %04d/%02d/%02d %02d:%02d:%02.0f H=%dkm (TECU)',...
          MjdToCal(td,ts+time(n)),hgts(1)))
    
    hold on, load('topo'), topo=topo([1,1:end,end],[181:end,1:181]);
    contour(-180:180,-90.5:90.5,topo,[0,0],'w')
    axis equal, axis([-180,180,-85,85])
    set(gca,'xtick',-180:30:180,'ytick',-90:30:90)
    xlabel('Longitude (deg)'), ylabel('Latitude (deg)')
end
