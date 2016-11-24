% “d—£‘w’x‰„ƒ‚ƒfƒ‹”äŠr
function sample7

dirs.pos='d:\gps\pos';
dirs.nav='d:\gps\nav';
dirs.ion='d:\gps\ion';

td=CalToMjd([2003,12,1]);
time=0:30:86400*30;
rcv='USUD';

posr=ReadPos(td,time(1),rcv,dirs.pos)';
gpos=EcefToGeod(posr);

azel=[0,90];
ion=repmat(nan,length(time),2);

for n=1:length(time)
    if mod(time(n),86400)==0
        [nav,inav,ionprm]=ReadNav(td,time(n),'','',dirs.nav);
    end
    ion(n,1)=ion_klob(td+time(n)/86400,azel*pi/180,gpos,ionprm);
    ion(n,2)=ion_tec(td,time(n),azel*pi/180,gpos,dirs.ion);
end

figure, plot(time/3600,ion), grid on
legend({'BROADCAST','IGS TEC MAP'})
xlabel('Time (H)'), ylabel('L1 ion-delay (m)')

dt1=MjdToCal(td+time(1)/86400); dt2=MjdToCal(td+time(end)/86400);
ti=sprintf('%04d/%02d/%02d %02d:%02d - %02d/%02d %02d:%02d',dt1(1:5),dt2(2:5));
title(sprintf('IONOSPHERE MODELS : %s : %s AZ=%.0f EL=%.0f',ti,rcv,azel(1),azel(2)))
