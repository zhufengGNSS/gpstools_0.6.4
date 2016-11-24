function [tt,zz,arc,i,slip,outl,loge,logs]=editobs(td,t,z,azel,sat,rcv,prm,f1,f2)
%-------------------------------------------------------------------------------
% [system] : GpsTools
% [module] : edit observation data
% [func]   : detect cycle-slip/outlier, repair slip and separate arcs
% [argin]  : td   = observation time (mjd)
%            t    = observation time vector (sec)
%            z    = observation data [L1,L2,P1,P2;...]
%            azel = azimath/elevation angle (rad)
%            dclk = clock jump offset (m)
%            sat  = satellite
%            rcv  = receiver
%           (prm) = edit parameters
%                prm.gapmax : max arc gap (sec)
%                prm.arcmin : min arc length (sec)
%                prm.sigmw  : noise of WM (m)
%                prm.sigif  : noise of ion-free phase (m)
%                prm.dgmax  : cycle-slip threshold of LG (m)
%                prm.dmp1   : cycle-slip threshold of MP1 (m)
%                prm.dmp2   : cycle-slip threshold of MP2 (m)
%                prm.dwl    : cycle-slip threshold of WL (cycle)
%                prm.outl   : outlier threshold (sigma)
%                prm.elwe   : threshold elevation weighting (1:on,0:off)
%                prm.wind   : slip detection window (points)
%                    prm.wind(1) : moving avarage window width
%                    prm.wind(2) : avaraging window width
%                prm.reps   : repair cycle-slip flag (1:on,0:off)
%                prm.npnt   : no of obs points for fitting
%                prm.nmax   : degree of poly. for fitting
%            (f1,f2) = L1,L2 frequency (Hz)
% [argout] : tt   = clean observation time vector (sec)
%            zz   = clean observation data
%            arc  = arc start/end index [arcstart1,arcend1;...]
%            i    = clean observation data index (zz=z(i,:))
%            slip = cycle-slip position [t,type;...]
%                       type 1-5:LG/MW/MP1/MP2/LC, 9:repair
%            outl = outlier position (sec)
%            loge,logs = edit slip/repair slip log
% [note]   :
% [version]: $Revision: 12 $ $Date: 2008-11-25 10:02:15 +0900 (火, 25 11 2008) $
%            Copyright(c) 2004-2006 by T.Takasu, all rights reserved
% [history]: 04/11/18  0.1  new
%            04/11/20  0.2  add argout nslip
%            05/06/24  0.3  add backward moving average check
%            05/08/14  0.4  add argin td,sat,rcv,edit parameter prm.elwe
%                           add argout slip,loge,logs
%            05/12/12  0.5  fix bug on detecting slip by moving avarage test
%            06/06/03  0.6  modify geometry-free screening
%-------------------------------------------------------------------------------
if nargin<4, prm=[]; end
if nargin<5, f1=1.57542E9; end
if nargin<6, f2=1.22760E9; end
if isempty(prm)
    prm.gapmax=300; prm.arcmin=600; prm.sigmw=0.5; prm.sigif=1.5;
    prm.dgmax=0.03; prm.dmp1=0.5; prm.dmp2=1; prm.dwl=0.5; prm.outl=4;
    prm.elwe=0; prm.wind=[32,16]; prm.reps=0; prm.npnt=10; prm.nmax=3;
end
C=299792458; lam1=C/f1; lam2=C/f2;
arc=[]; t=t(:); slip=[]; outl=[]; loge={}; logs={}; rep=[];

% separate arcs
[tt,zz,arc0,i]=SeparateArc(t,z,prm.gapmax,prm.arcmin);
if isempty(arc0), return, end

% detect slip by LG difference
[arc1,s1,sig1]=ScreenLG(tt,zz,arc0,azel(i,2),prm.dgmax,lam1,lam2,prm.elwe);

% detect slip by MW difference
[arc2,s2,sig2]=ScreenMW(tt,zz,arc1,azel(i,2),prm.sigmw,prm.dwl,lam1,lam2,prm.wind,prm.elwe);

% detect slip by MP1 differenc
[arc3,s3]=ScreenMP1(tt,zz,arc2,azel(i,2),prm.dmp1,lam1,lam2,prm.wind,prm.elwe);

% detect slip by MP2 difference
[arc4,s4]=ScreenMP2(tt,zz,arc3,azel(i,2),prm.dmp2,lam1,lam2,prm.wind,prm.elwe);

% detect slip by LC-PC difference
[arc5,s5]=ScreenLC_PC(tt,zz,arc4,azel(i,2),prm.sigif,lam1,lam2,prm.wind,prm.elwe);

slip=sortrows([tt(s1),repmat(1,length(s1),1);tt(s2),repmat(2,length(s2),1);...
               tt(s3),repmat(3,length(s3),1);tt(s4),repmat(4,length(s4),1);...
               tt(s5),repmat(5,length(s5),1)]);
% repair slips
if prm.reps
    [zz,arc5,rep,logs]=RepairSlip(td,tt,zz,arc0,arc5,azel(i,2),prm.npnt,...
                                  prm.nmax,sat,rcv,lam1,lam2);
    [ts,j]=intersect(slip(:,1),tt(rep)); slip(j,2)=9;
end
% delete outliers
if prm.outl>0
    [zz,outl]=RemoveOutlier(tt,zz,arc5,azel(i,2),prm.outl,lam1,lam2);
end
% delete arcs shorter than min arc length
if ~isempty(arc5)
    tz=round(tt); arc=arc5(tz(arc5(:,2))-tz(arc5(:,1))>=prm.arcmin,:);
end
ns=[length(s1),length(s2),length(s3),length(s4),length(s5)];
loge={sprintf('%-7s %-7s:%6d %3d %3d %3d %3d %3d %3d %3d  %3d  %3d %3d %6.3f %6.3f',...
              sat,rcv,size(zz,1),ns,sum(ns),length(rep),length(outl),size(arc0,1),...
              size(arc,1),sig1*100,sig2)};

% separate arcs ----------------------------------------------------------------
function [t,z,arc,i]=SeparateArc(t,z,gapmax,arcmin)
arc=[]; i=find(all(~isnan(z),2)); if isempty(i), return, end
t=t(i); z=z(i,:); tz=round(t); j=find(diff(tz)>gapmax);
arc=[1,j'+1;j',length(t)]';
if ~isempty(arc), arc=arc(tz(arc(:,2))-tz(arc(:,1))>=arcmin,:); end

% detect cycle-slip by LG double-difference ------------------------------------
function [arc,slip,sig]=ScreenLG(t,z,arc,el,dgmax,lam1,lam2,elwe)
narc=[]; slip=[]; sig=0; n=0; if dgmax<=0, return, end
for a=arc'
    i=a(1):a(2);
    if elwe, thres=dgmax./sin(el(i(2:end-1))); else thres=dgmax; end
    j=find(abs(diff(diff(LG(z(i,:),lam1,lam2))))>thres);
    k=unique([j,j+1]); k(k>length(i)-1)=[];
    narc=[narc;[i(1),i(k+1);i(k),i(end)]']; slip=[slip,i(k+1)];
end
arc=narc;
for a=narc'
    i=a(1):a(2); m=length(i);
    if m>2, s=diff(diff(LG(z(i,:),lam1,lam2))); sig=sig+s'*s; n=n+m-2; end
end
if n>0, sig=sqrt(sig/n)/2; end

% detect slip by MW difference -------------------------------------------------
function [arc,slip,sig]=ScreenMW(t,z,arc,el,sigmw,dwl,lam1,lam2,wind,elwe)
lam4=1/(1/lam1-1/lam2); lam5=1/(1/lam1+1/lam2);
narc=[]; slip=[]; sig=0; n=0; if dwl<=0, return, end
for a=arc'
    i=a(1):a(2);
    j=find(abs(diff(diff(LG(z(i,:),lam1,lam2))))>0.05)';
    k=unique([j+1,j+2]); k(k>length(i)-1)=[];
    j=DetectSlip(t(i),MW(z(i,:),lam1,lam2,lam4,lam5),k,el(i),lam4*dwl,sigmw,...
                 wind,elwe);
    narc=[narc;[i(1),i(j);i(j-1),i(end)]']; slip=[slip,i(j)];
end
arc=narc;
for a=arc'
    i=a(1):a(2); m=length(i);
    if m>2, s=diff(diff(MW(z(i,:),lam1,lam2,lam4,lam5))); sig=sig+s'*s; n=n+m-2; end
end
if n>0, sig=sqrt(sig/n)/2; end

% detect slip by MP1 difference ------------------------------------------------
function [arc,slip]=ScreenMP1(t,z,arc,el,dmp,lam1,lam2,wind,elwe)
narc=[]; slip=[]; if dmp<=0, return, end
for a=arc'
    i=a(1):a(2);
    j=DetectSlip(t(i),MP1(z(i,:),lam1,lam2),[],el(i),dmp,dmp*2,wind,elwe);
    narc=[narc;[i(1),i(j);i(j-1),i(end)]']; slip=[slip,i(j)];
end
arc=narc;

% detect slip by MP2 difference ------------------------------------------------
function [arc,slip]=ScreenMP2(t,z,arc,el,dmp,lam1,lam2,wind,elwe)
narc=[]; slip=[]; if dmp<=0, return, end
for a=arc'
    i=a(1):a(2);
    j=DetectSlip(t(i),MP2(z(i,:),lam1,lam2),[],el(i),dmp,dmp*2,wind,elwe);
    narc=[narc;[i(1),i(j);i(j-1),i(end)]']; slip=[slip,i(j)];
end
arc=narc;

% detect slip by LC-PC differences ---------------------------------------------
function [arc,slip]=ScreenLC_PC(t,z,arc,el,sigif,lam1,lam2,wind,elwe)
narc=[]; slip=[]; if sigif<=0, return, end
for a=arc'
    i=a(1):a(2);
    j=DetectSlip(t(i),LC_PC(z(i,:),lam1,lam2),[],el(i),sigif*4,sigif,wind,elwe);
    narc=[narc;[i(1),i(j);i(j-1),i(end)]']; slip=[slip,i(j)];
end
arc=narc;

% detect slip by moving average test -------------------------------------------
function slip=DetectSlip(t,z,slip,el,thres,sig,wind,elwe)
i=1; ss=[]; slip=[]; if any(wind<=0), return, end
while ~isempty(i) % forward moving average test
    if elwe, s=sig/sin(el(i)); else s=sig; end
    [j,sslip]=MovingAveTest(z(i:end),wind(1),s);
    i=i+j-1; slip=[slip,i]; if sslip>10, ss=[ss,i]; end
end
i=length(z); 
while ~isempty(i) % backward moving average test
    if elwe, s=sig/sin(el(i)); else s=sig; end
    [j,sslip]=MovingAveTest(z(i:-1:1),wind(1),s);
    i=i-j+1; slip=[slip,i+1]; if sslip>10, ss=[ss,i+1]; end
end
slip=unique(slip); k=[];
for n=1:length(slip) % identify slip if difference of avarages over threshold
    i=slip(n)-wind(2); j=slip(n)+wind(2)-1;
    i=max([i,ss(ss<slip(n)),1]); 
    j=min([j,ss(ss>slip(n))-1,length(z)]); 
    i=i:slip(n)-1; j=slip(n):j;
    np=length(i); nf=length(j);
    if abs(sum(z(i))/np-sum(z(j))/nf)>thres*(1+1/sqrt(min(np,nf))), k=[k,n]; end
end
slip=slip(k);

% repair slip ------------------------------------------------------------------
function [z,arc,rep,log]=RepairSlip(td,t,z,arc0,arc,el,npnt,nmax,sat,rcv,lam1,lam2)
lam4=1/(1/lam1-1/lam2); lam5=1/(1/lam1+1/lam2);

% generate LG-slip table
table=[];
for n1=-16:16, for n2=-16:16, table=[table;n1,n2,n1-n2,lam1*n1-lam2*n2]; end, end
table=sortrows(table);

narc=[]; rep=[]; log={};
for a=arc0'
    narc=[narc;a(1:2)'];
    for i=min(find(a(1)<=arc(:,1))):max(find(a(2)>=arc(:,2)))-1
        
        slip=arc(i+1,1); fix=0;
        j=arc(i,1):arc(i+1,2); zi=z(j,:); ti=t(j); si=slip-j(1)+1; 
        
        % fix guess slip of L1/L2
        n1=round(EstSlipF(ti,zi(:,1),si,2,1));
        n2=round(EstSlipF(ti,zi(:,2),si,2,1));
        zi(si:end,1)=zi(si:end,1)-n1;
        zi(si:end,2)=zi(si:end,2)-n2;
        
        % estimate slip amount by MW and LG
        mw=EstSlipM(ti,MW(zi,lam1,lam2,lam4,lam5),si,npnt)/lam4;
        gf=EstSlipF(ti,LG(zi,lam1,lam2),si,npnt,nmax);
        
        if ~isempty(mw)&~isempty(gf)&abs(mw-round(mw))<0.3
            
            % search LG-slip table
            tbl=table(table(:,3)==round(mw),:);
            [s,j]=sort(abs(tbl(:,4)-gf));
            
            if length(s)>=2&s(1)*3<s(2) % ratio-test
                n1=n1+tbl(j(1),1); n2=n2+tbl(j(1),2); fix=1;
                z(slip:a(2),1)=z(slip:a(2),1)-n1; % fix slip
                z(slip:a(2),2)=z(slip:a(2),2)-n2;
                msg=sprintf('%-7s %-7s: %s %9d %9d %6.2f %5.1f %4.1f %4.1f',sat,...
                            rcv,tstr(td,t(slip)),n1,n2,mw,gf*100,s(1:2)*100);
                rep=[rep,slip]; log={log{:},msg};
            end
        end
        if ~fix
            narc(end+1,:)=[arc(i+1,1),narc(end,2)];
            narc(end-1,:)=[narc(end-1,1),arc(i,2)];
        end
    end
end
arc=narc;

% delete outliers --------------------------------------------------------------
function [z,outl]=RemoveOutlier(t,z,arc,el,thres,lam1,lam2)
lam4=1/(1/lam1-1/lam2); lam5=1/(1/lam1+1/lam2); outl=[];
for a=arc'
    i=a(1):a(2);
    z1=MW(z(i,:),lam1,lam2,lam4,lam5);
    z2=LC_PC(z(i,:),lam1,lam2);
    z1=(z1-mean(z1)).*sin(el(i));
    z2=(z2-mean(z2)).*sin(el(i));
    s1=std(z1)*thres;
    s2=std(z2)*thres;
    j=i(z1<-s1|s1<z1|z2<-s2|s2<z2);
    z(j,:)=nan; outl=[outl;t(j)];
end

% moving average test ----------------------------------------------------------
function [slip,sslip]=MovingAveTest(z,wind,sig)
slip=[]; sslip=0;
mave=zeros(size(z)); mvar=zeros(size(z)); mave(1)=z(1); mvar(1)=sig^2;
for n=2:length(z)
    sslip=abs(z(n)-mave(n-1))/sqrt(mvar(n-1)); if sslip>4, slip=n; break; end
    if n<=wind
        mave(n)=mave(n-1)+(z(n)-mave(n-1))/n;
        mvar(n)=(sum((z(1:n)-mave(n)).^2)+sig^2*(wind-n))/wind;
    else
        mave(n)=mave(n-1)+(z(n)-z(n-wind))/wind;
        mvar(n)=sum((z(n-wind+1:n)-mave(n)).^2)/wind;
    end
end

% estimate slip amount by difference of averages -------------------------------
function sslip=EstSlipM(t,z,slip,npnt)
if slip-1<npnt|length(z)-slip+1<npnt, sslip=[]; return, end
sslip=mean(z(slip:end))-mean(z(1:slip-1));

% estimate slip amount by polynominal fitting ----------------------------------
function sslip=EstSlipF(t,z,slip,npnt,nmax)
sslip=[]; i=max(1,slip-npnt):min(length(z),slip+npnt-1);
if length(i)<nmax+1, return, end
t=t(i); t=(2*t-(t(end)+t(1)))/(t(end)-t(1)); % normalize to (-1,1)
A=ones(length(i),nmax+2); A(1:slip-i(1),nmax+2)=0;
for n=nmax:-1:1, A(:,n)=A(:,n+1).*t; end
c=A\z(i); sslip=c(end);

% MW (Melbourne-Wubbena) combination -------------------------------------------
function zz=MW(z,lam1,lam2,lam4,lam5)
zz=lam4*(z(:,1)-z(:,2))-lam5*(z(:,3)/lam1+z(:,4)/lam2);

% LG (geometry-free phase) combination -----------------------------------------
function zz=LG(z,lam1,lam2)
zz=lam1*z(:,1)-lam2*z(:,2);

% MP1 (L1 multipath) combination -----------------------------------------------
function zz=MP1(z,lam1,lam2)
c2=lam1^2/(lam2^2-lam1^2);
zz=z(:,3)-lam1*z(:,1)-2*c2*(lam1*z(:,1)-lam2*z(:,2));

% MP2 (L2 multipath) combination -----------------------------------------------
function zz=MP2(z,lam1,lam2)
c1=lam2^2/(lam2^2-lam1^2);
zz=z(:,4)-lam2*z(:,2)-2*c1*(lam1*z(:,1)-lam2*z(:,2));

% LC-PC (ion-free phase-code) combination --------------------------------------
function zz=LC_PC(z,lam1,lam2)
c1=lam2^2/(lam2^2-lam1^2); c2=lam1^2/(lam2^2-lam1^2);
zz=c1*lam1*z(:,1)-c2*lam2*z(:,2)-c1*z(:,3)+c2*z(:,4);

% time string ------------------------------------------------------------------
function s=tstr(td,t), s=sprintf('%04d/%02d/%02d %02d:%02d:%02.0f',mjdtocal(td,t));
