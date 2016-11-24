function t_sitedisp % ŠÏ‘ª‹ÇˆÊ’u•ÏˆÊ’P‘ÌŽŽŒ±

td=caltomjd([2003,12,1]);
time=0:1800:86400;
%rcvs={'AMC2','ALBH','SANT','HRAO','PERT','TSKB'};
%rcvs={'AMC2'};
rcvs={'TSKB'};
dirs.pos='d:\gps\pos';
dirs.erp='d:\gps\erp';
utc_tai=-32;
opt=[1,0,0,1];

for n=1:length(rcvs)
    disp('SITE DISPLACEMENT(m)')
    disp(sprintf('station    = %s',rcvs{n}))
    disp(sprintf('start time = %04d/%02d/%02d %02d:%02d:%02d',mjdtocal(td)))
    posr=readpos(td,time(1),rcvs{n},dirs.pos)';
    gpos=eceftogeod(posr);
    [epos,E]=geodtoecef(gpos);
    [odisp,ophas]=readoload(rcvs{n},dirs.pos,'igs');
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
    figure, plot(time/3600,dpf), legend({'E(m)','N(m)','R(m)'})
    grid on, xlabel('time(H)')
    title(['SITE DISPLACEMENT : ',rcvs{n}])
    
    dt=mjdtocal(td);
    load(fullfile(dirs.pos,sprintf('sdsp_%s_%04d%02d%02d%02d%02d.mat',rcvs{n},dt(1:5))))
    figure, plot(time/3600,dposf),grid on,title('GOTIC2')
    figure, plot(time/3600,dposf-dpf),grid on,title('GOTIC2-SITEDISP')
end
