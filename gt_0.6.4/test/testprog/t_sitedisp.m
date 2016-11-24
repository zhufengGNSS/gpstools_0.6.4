function t_sitedisp % ŠÏ‘ª‹ÇˆÊ’u•ÏˆÊ’P‘ÌŽŽŒ±

td=caltomjd([2007,1,1]);
time=0:1800:86400*31;
%rcvs={'AMC2','ALBH','SANT','HRAO','PERT','TSKB'};
%rcvs={'AMC2'};
rcvs={'TSKB'};
dirs.pos='k:\gps\pos';
dirs.erp='k:\gps\erp';
utc_tai=-33;
opt=[1,1,1,0];
file='../../data/oload_igs.blq';

for n=1:length(rcvs)
    disp('SITE DISPLACEMENT(m)')
    disp(sprintf('station    = %s',rcvs{n}))
    disp(sprintf('start time = %04d/%02d/%02d %02d:%02d:%02d',mjdtocal(td)))
    posr=readpos(td,time(1),rcvs{n},dirs.pos,'approx')';
    gpos=eceftogeod(posr);
    [epos,E]=geodtoecef(gpos);
    [odisp,ophas]=readoload(rcvs{n},file);
    dp=zeros(length(time),3);
    dpf=zeros(length(time),3);
    for m=1:length(time)
        tutc=td+(time(m)+19+utc_tai)/86400;
        [rsun,rmoon]=sunmoonpos(tutc);
        erp=readerp(tutc,dirs.erp);
        [U,P,N,gmst]=ecsftoecef(tutc,erp,utc_tai);
        dp(m,:)=sitedisp(tutc,posr,U*rsun,U*rmoon,odisp,ophas,gmst,erp,opt)';
        dpf(m,:)=(E*dp(m,:)')';
    end
    label={'E-W (cm)','N-S (cm)','U-D (cm)'};
    figure('color','w')
    for m=1:3
        ggt('subplotv','',m,3,'margin',[.08,.005,.08,.01],'fontsize',12);
        plot(time/3600,dpf(:,m)*100,'k','linewidth',1.5)
        xlim([time(1),time(end)]/3600);
        ylim([-40,40]);
        ylabel(label{m});
        if m==3, xlabel('Time (H)'); else set(gca,'xticklabel',[]); end
    end
end
