function t_pos

td=caltomjd([2007,1,1]);
load('d:\test\posf_TSKB_2007010100.mat')
dirs.pos='k:\gps\pos';
dirs.erp='k:\gps\erp';

pos=readpos(td,0,'TSKB',dirs.pos,'igssnx')';
gpos=eceftogeod(pos);
[e,E]=geodtoecef(gpos);

for m=1:length(time)
    err(m,:)=(E*(data(m,1:3)'-pos))';
    sig(m,:)=sqrt(diag(E*diag(covs(m,1:3))*E')');
end
label={'E-W (cm)','N-S (cm)','U-D (cm)'};
figure('color','w')
for m=1:3
    ggt('subplotv','',m,3,'margin',[.08,.005,.08,.01],'fontsize',12,...
        'xtick',0:2:24,'ytick',-10:2:10);
    plot(time/3600,err(:,m)*100,'k.-','markersize',8)
    plot(time/3600,sig(:,m)*100,'k--','linewidth',0.3)
    plot(time/3600,-sig(:,m)*100,'k--','linewidth',0.3)
    xlim([0,24]);
    ylim([-5,5]);
    ylabel(label{m});
    if m==3, xlabel('Time (H)'); else set(gca,'xticklabel',[]); end
end
