function gpsmpc(prm)
%-------------------------------------------------------------------------------
% [system] : GpsTools
% [module] : code multipath profiler
% [func]   : estimate code multipath profile
% [argin]  :(prm) = parameters struct
% [argout] : none
% [note]   : output file
%            mpc_{rcv}_YYYYMMDDHH.mat : code multipath profile
% [version]: $Revision: 8 $
%            Copyright(c) 2004-2006 by T.Takasu, all rights reserved
% [history]: 04/10/25  0.1  new
%-------------------------------------------------------------------------------
global p_, if nargin<1, p_=prm_gpsest; else p_=prm; end

for n=1:length(p_.rcvs)
    EstMpc(p_.td,p_.time,p_.sats,p_.rcvs{n},p_.dirs,p_.mpc.nmax,...
           p_.mpc.narcmax,p_.mpc.elmin);
end
clear p_

% estimate multipath -----------------------------------------------------------
function EstMpc(td,time,sats,rcv,dirs,nmax,narcmax,elmin)

[z,iz]=ReadObs(td,time,sats,rcv,dirs.obs); if isempty(iz), return, end
posr=ReadPos(td,time(1),rcv,dirs.pos)';

nx=(nmax+1)^2-1;
AA=zeros(nx+narcmax); Ay1=zeros(nx+narcmax,1); Ay2=Ay1; narc=0;

for n=1:length(sats)
    i=find(iz(:,2)==n); izz=iz(i,:); zz=z(i,:);
    if ~isempty(i)
        nav=readnav(td,time,'',dirs.nav);
        azel=zeros(length(i),2);
        for m=1:size(izz,1)
            azel(m,:)=SatAzEl(NavToState(td,izz(m,1),nav),posr);
        end
        [t,zz,arc,i]=EditObs(izz,zz,azel); azel=azel(i,:); azels{n}=[];
        for a=arc'
            narc=narc+1; i=a(1):a(2);
            [mp1,mp2,j]=MultiPath(zz(i,:));
            disp(sprintf('%s-%s : arc=%7d-%7dsec std(mp1)=%6.3fm std(mp2)=%.3fm',...
                 sats{n},rcv,t(a(1)),t(a(2)),std(mp1),std(mp2)))
            sig=1./sin(azel(i(j),2));
            A=sparse(size(mp1,1),narcmax); A(:,narc)=1;
            A=[DesignMat(azel(i(j),:),nmax),A];
            A=A./repmat(sig,1,size(A,2));
            AA=AA+A'*A;
            Ay1=Ay1+A'*(mp1./sig);
            Ay2=Ay2+A'*(mp2./sig);
            azels{n}=[azels{n};azel(i(j),:)];
        end
    end
end
x1=repmat(nan,nx+narcmax,1); x2=x1;
i=find(any(AA));
x1(i)=AA(i,i)\Ay1(i);
x2(i)=AA(i,i)\Ay2(i);
[mpcc1,mpcs1]=DeserialSphFunc(x1(1:nx),nmax);
[mpcc2,mpcs2]=DeserialSphFunc(x2(1:nx),nmax);
mpcc1(1,1)=-MpOffset(mpcc1,mpcs1,azels,elmin);
mpcc2(1,1)=-MpOffset(mpcc2,mpcs2,azels,elmin);

time=time'; sats=sats'; dt=MjdToCal(td);
file=fullfile(dirs.est,sprintf('mpcs_%s_%04d%02d%02d.mat',rcv,dt(1:3)));
save(file,'td','time','sats','rcv','mpcc1','mpcs1','mpcc2','mpcs2','azels','nmax')
disp(['save : ',file])

% generate design matrix -------------------------------------------------------
function A=DesignMat(azel,nmax)
A=zeros(size(azel,1),(nmax+1)^2-1);
for n=1:size(azel,1)
    [fc,fs]=SphFunc(azel(n,1),azel(n,2),nmax);
    A(n,:)=SerialSphFunc(fc,fs,nmax)';
end

% calculate code multipath by phase obs. ---------------------------------------
function [mp1,mp2,i]=MultiPath(z)
C=299792458; f1=1.57542E9; f2=1.22760E9; lam1=C/f1; lam2=C/f2;
mp1=repmat(nan,size(z,1),1); mp2=mp1;
for n=1:size(z,1)
    cp1=lam1*z(n,1);
    cp2=lam2*z(n,2);
    ion=-(cp1-cp2)/(1-f1^2/f2^2);
    mp1(n)=z(n,3)-cp1-2*ion;
    mp2(n)=z(n,4)-cp2-2*f1^2/f2^2*ion;
end
i=find(~isnan(mp1)&~isnan(mp1));
mp1=mp1(i)-mean(mp1(i));
mp2=mp2(i)-mean(mp2(i));

% calculate offset as avarage of multipaths ------------------------------------
function off=MpOffset(mpcc,mpcs,azels,elmin)
off=0; nmp=0;
for n=1:length(azels)
    for m=find(azels{n}(:,2)>=elmin)'
        off=off+RcvMpc(azels{n}(m,:),mpcc,mpcs); nmp=nmp+1;
    end
end
if nmp>0, off=off/nmp; end

% serialize spherical harmonic coefficients ------------------------------------
function f=SerialSphFunc(c,s,nmax)
f=zeros((nmax+1)^2-1,1); i=1;
for n=1:nmax, for m=0:n
    if m==0, f(i)=c(n+1,1); i=i+1; else f(i:i+1)=[c(n+1,m+1);s(n+1,m+1)]; i=i+2; end
end, end

% deserialize spherical harmonic coefficients ----------------------------------
function [c,s]=DeserialSphFunc(f,nmax)
c=zeros(nmax+1,nmax+1); s=c; i=1;
for n=1:nmax, for m=0:n
    if m==0, c(n+1,1)=f(i); i=i+1; else c(n+1,m+1)=f(i); s(n+1,m+1)=f(i+1); i=i+2; end
end, end
