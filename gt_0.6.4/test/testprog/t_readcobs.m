td=caltomjd([2003,12,1]);
%time=0:86400;
time=3600*16:3600*24;
obsdir='obs\rinexh';
obssrc='rinexh';
navdir='nav';
posdir='pos';
rcv='MATE';

for n=1:32, sats{n}=sprintf('GPS%02d',n); end
%[z,iz]=readobs(td,time,sats,rcv,obsdir,obssrc);
%i=find(iz(:,2)==3);
%figure, plot(iz(i,1)/3600,z(i,1:2));

obs=readcleanobs(td,time,sats,rcv,obsdir,obssrc,navdir,posdir);

