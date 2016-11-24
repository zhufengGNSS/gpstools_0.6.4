function [x,ix,P,sig,prn]=gt_initstate(prm)
%-------------------------------------------------------------------------------
% [system] : GpsTools
% [module] : initialize state variables
% [func]   : initialize state variables
% [argin]  : prm = processing parameters struct (see prm_gpsest_def.m)
% [argout] : x   = state variables vector
%            ix  = state variables index
%                ix.sato{s}   = satellite s orbit
%                ix.sats{s}   = satellite s orbit parameters
%                ix.satc{s}   = satellite s clock
%                ix.sata{s}   = satellite s all states
%                ix.rcvc{r}   = station r reciever clock
%                ix.rcvz{r}   = station r tropospheric parameters
%                ix.rcvp{r}   = station r position
%                ix.arcn{s,r} = satellite s - station r phase bias
%                ix.rcva{r}   = station r all states (include phase bias)
%                ix.erp/eco   = earth rotation parameters/geocenter offset
%            P   = state variables covariance matrix
%            sig = apriori state std. devs. vector
%            prn = process noise std. devs. vector
% [note]   :
% [version]: $Revision: 5 $ $Date: 06/07/08 1:16 $
%            Copyright(c) 2004-2006 by T.Takasu, all rights reserved
% [history]: 06/03/26  0.1  separated from gpsestd.m
%-------------------------------------------------------------------------------
[ix,nx]=initindex(prm);
x=repmat(nan,nx,1);
P=zeros(nx);
sig=initsig(ix,nx,prm);
prn=initprn(ix,nx,prm);

% initialize state index -------------------------------------------------------
function [ix,n]=initindex(prm)
ix=[]; n=0; ix.wn=[];
for m=1:length(prm.sats)
    ix.sato{m}=[]; ix.sats{m}=[]; ix.satc{m}=[];
    if prm.est.sato(m)==1, ix.sato{m}=n+(1:6); n=n+6; end
    if prm.est.sats(m)==1, ix.sats{m}=n+(1:6); n=n+6; end
    if prm.est.satc(m)==1
        switch prm.sclkmodel
        case 0, ix.satc{m}=n+1; n=n+1;
        case 1, ix.satc{m}=n+(1:2); n=n+2; end
    end
    ix.sata{m}=[ix.sato{m},ix.sats{m},ix.satc{m}];
end
for m=1:length(prm.rcvs)
    ix.rcvc{m}=[]; ix.rcvz{m}=[]; ix.rcvg{m}=[]; ix.rcvp{m}=[];
    if prm.est.rcvc(m)==1
        switch prm.rclkmodel
        case 0, ix.rcvc{m}=n+1; n=n+1;
        case 1, ix.rcvc{m}=n+(1:2); n=n+2; end
    end
    if prm.est.rcvz(m)==1
        switch prm.zpdmodel
        case 0, ix.rcvz{m}=n+(1:1); n=n+1;
        case 1, ix.rcvz{m}=n+(1:2); n=n+2; end
        switch prm.trgmodel
        case 1, ix.rcvg{m}=n+(1:2); n=n+2;
        case 2, ix.rcvg{m}=n+(1:5); n=n+5; end
    end
    if prm.est.rcvp(m)==1
        switch prm.rposmodel
        case {0,1}, ix.rcvp{m}=n+(1:3); n=n+3;
        case 2,     ix.rcvp{m}=n+(1:6); n=n+6; end
    end
end
for m=1:length(prm.sats)
    for k=1:length(prm.rcvs)
        ix.arcn{m,k}=[];
        if prm.est.rcvb(k)==0, ix.arcn{m,k}=n+1; n=n+1; end
    end
end
for m=1:length(prm.rcvs)
    ix.rcva{m}=[ix.rcvc{m},ix.rcvz{m},ix.rcvg{m},ix.rcvp{m},ix.arcn{:,m}];
end
if prm.est.erp==1, ix.erp=n+(1:3); n=n+3; else ix.erp=[]; end
if prm.est.eco==1, ix.eco=n+(1:4); n=n+4; else ix.eco=[]; end

% initialize apriori state std. dev. -------------------------------------------
function sig=initsig(ix,nx,prm)
sig=zeros(1,nx);
for n=1:length(prm.sats)
    i=ix.sato{n}; sig(i)=prm.sig.sato(n,1:length(i));
    i=ix.sats{n}; sig(i)=prm.sig.sats(n,1:length(i));
    i=ix.satc{n}; sig(i)=prm.sig.satc(n,1:length(i));
end
for n=1:length(prm.rcvs)
    i=ix.rcvc{n}; sig(i)=prm.sig.rcvc(n,1:length(i));
    i=ix.rcvz{n}; sig(i)=prm.sig.rcvz(n,1:length(i));
    i=ix.rcvg{n}; sig(i)=prm.sig.rcvg(n,1:length(i));
    i=ix.rcvp{n}; sig(i)=prm.sig.rcvp(n,1:length(i));
end
for n=1:length(prm.sats),
    for m=1:length(prm.rcvs)
        i=ix.arcn{n,m}; sig(i)=prm.sig.arcn(m,1:length(i));
    end
end
if prm.est.erp==1, i=ix.erp; sig(i)=prm.sig.erp; end
if prm.est.eco==1, i=ix.eco; sig(i)=prm.sig.eco; end

% initialize process noise std. dev. -------------------------------------------
function prn=initprn(ix,nx,prm)
prn=zeros(1,nx);
for n=1:length(prm.sats)
    i=ix.sato{n}; prn(i)=prm.prn.sato(n,1:length(i));
    i=ix.sats{n}; prn(i)=prm.prn.sats(n,1:length(i));
    i=ix.satc{n}; prn(i)=prm.prn.satc(n,1:length(i));
end
for n=1:length(prm.rcvs)
    i=ix.rcvc{n}; prn(i)=prm.prn.rcvc(n,1:length(i));
    i=ix.rcvz{n}; prn(i)=prm.prn.rcvz(n,1:length(i));
    i=ix.rcvg{n}; prn(i)=prm.prn.rcvg(n,1:length(i));
    i=ix.rcvp{n}; prn(i)=prm.prn.rcvp(n,1:length(i));
end
for n=1:length(prm.sats),
    for m=1:length(prm.rcvs)
        i=ix.arcn{n,m}; prn(i)=prm.prn.arcn(m,1:length(i));
    end
end
if prm.est.erp==1, i=ix.erp; prn(i)=prm.prn.erp; end
if prm.est.eco==1, i=ix.eco; prn(i)=prm.prn.eco; end
