function gpsmpp(td,time,rcvs,dirs,outdir,fb,tunit,nmax)
%-------------------------------------------------------------------------------
% [system] : GpsTools
% [module] : phase multipath profiler
% [func]   : estimate phase multi path profile
% [argin]  : td     = start time (mjd-gpst)
%            time   = time vector
%            rcvs   = stations
%           (dirs)  = directory
%               dirs.est = estimation data direcotry
%               dirs.nav = navigation data direcotry
%           (outdir) = outputdirectory
%           (fb)    = estimation direction ('f':forward,'b':backward)
%           (tunit) = processing unit time
%           (nmax)  = degree of spheric harmonic function
% [argout] : none
% [note]   : output file
%            mpc_{rcv}_YYYYMMDDHH.mat : phase multipath profile
% [version]: $Revision: 12 $
%            Copyright(c) 2004-2006 by T.Takasu, all rights reserved
% [history]: 04/11/19  0.1  new
%            05/06/10  0.2  use postfit-residuals
%-------------------------------------------------------------------------------
if nargin<4, dirs.est=''; dirs.nav=''; dirs.pos=''; end
if nargin<5, outdir=''; end
if nargin<6, fb='f'; end
if nargin<7, tunit=24; end
if nargin<8, nmax=8; end

for n=1:length(rcvs)
    [t,ress,azel,index,sats]=ReadRess(td,time,rcvs(n),dirs,fb,tunit);
    if ~isempty(t)
        EstMpp(td,t,ress,azel,index(:,1:2),sats,rcvs{n},outdir,nmax);
    end
end

% estimate multipath by residuals ----------------------------------------------
function EstMpp(td,t,ress,azel,index,sats,rcv,outdir,nmax)
nx=(nmax+1)^2-1;
AA=zeros(nx); Ay=zeros(nx,1);
for i=find(~isnan(ress))'
    A=DesignMat(azel(i,:),nmax);
    AA=AA+A'*A;
    Ay=Ay+A'*ress(i);
end
[mpcc,mpcs]=DeserialSphFunc(AA\Ay,nmax);
mpcc=[zeros(1,size(mpcc,2));mpcc];
mpcs=[zeros(1,size(mpcc,2));mpcs];

time=unique(t); sats=sats'; dt=MjdToCal(td); azels={};
for n=1:length(sats), azels{n}=azel(index(:,1)==n,:); end
file=gfilepath(outdir,sprintf('mpps_%s_%04d%02d%02d.mat',rcv,dt(1:3)),dt);
save(file,'td','time','sats','rcv','mpcc','mpcs','azels','nmax')
disp(['save : ',file])

% design matrix ----------------------------------------------------------------
function A=DesignMat(azel,nmax)
A=zeros(size(azel,1),(nmax+1)^2-1);
for n=1:size(azel,1)
    [fc,fs]=SphFunc(azel(n,1),azel(n,2),nmax);
    A(n,:)=SerialSphFunc(fc(2:end,:),fs(2:end,:),nmax)';
end

% serialize spheric harmonic funcs ---------------------------------------------
function f=SerialSphFunc(c,s,nmax)
f=zeros((nmax+1)^2-1,1); i=1;
for n=1:nmax, for m=0:n
    if m==0, f(i)=c(n,1); i=i+1; else f(i:i+1)=[c(n,m+1);s(n,m+1)]; i=i+2; end
end, end

% deserialize spheric harmonic funcs -------------------------------------------
function [c,s]=DeserialSphFunc(f,nmax)
c=zeros(nmax,nmax+1); s=c; i=1;
for n=1:nmax, for m=0:n
    if m==0, c(n,1)=f(i); i=i+1; else c(n,m+1)=f(i); s(n,m+1)=f(i+1); i=i+2; end
end, end

% read postfit residuals ------------------------------------------------------
function [t,ress,azel,index,sats]=ReadRess(td,time,rcvs,dirs,fb,tunit);
t=[]; ress=[]; azel=[]; index=[]; sats={};
tu=tunit*3600;
for ts=(floor(time(1)/tu):floor(time(end)/tu))*tu
    [tn,indn,resn,outn,sats]=...
        readstats(td,time(ts<=time&time<ts+tu),dirs.est,['res',fb],tunit,rcvs);
    i=find(~outn&ts<=tn&tn<ts+tu);
    if ~isempty(i)
        t=[t;tn(i)]; ress=[ress;resn(i,2)]; index=[index;indn(i,:)];
        azel=[azel;SatDir(td,tn(i),indn(i,:),sats,rcvs,dirs)];
    end
end

% calculate satellite/station directions ---------------------------------------
function azels=SatDir(td,time,index,sats,rcvs,dirs)
[nav,inav]=readnav(td,time(1):86400:time(end-1),sats,'',dirs.nav);
posr=shiftdim(readpos(td,time(1),rcvs,'','approx'),1);
[t,it,jt]=unique(time);
poss=zeros(3,length(t),length(sats)); possun=zeros(3,length(t));
for n=1:length(sats)
    navs=nav(find(inav==n),:);
    for m=1:length(t)
        poss(:,m,n)=navtostate(td,t(m),navs);
    end
end
azels=zeros(length(time),2);
for n=1:length(time)
    i=index(n,1); j=index(n,3);
    azels(n,:)=satazel(poss(:,jt(n),i),posr(:,j));
end
