function [state,abort,msg,stats]=gt_readstate(td,time,state,prm,opt)
%-------------------------------------------------------------------------------
% [system] : GpsTools
% [module] : read a priori states
% [func]   : read a priori states
% [argin]  : td,time = day(mjd), time(sec)
%            state   = a priori states
%            prm     = processing parameters (see prm_gpsest_def.m)
%            (opt)   = option (1:satellite states,2:station states,0:both)
%                      (default:0)
% [argout] : state   = predetermined states
%                state(n).ephs(:,s) = time(n) prm.sats{s} pos/vel
%                state(n).ephv(:,s) = time(n) prm.sats{s} pos/vel variance
%                state(n).clks(:,s) = time(n) prm.sats{s} clock
%                state(n).clkv(:,s) = time(n) prm.sats{s} clock variance
%                state(n).clkr(:,r) = time(n) prm.rcvs{r} clock
%                state(n).zpdr(:,r) = time(n) prm.rcvs{r} tropospheric ztd
%                state(n).posr(:,r) = time(n) prm.rcvs{r} position
%                state(n).posv(:,r) = time(n) prm.rcvs{r} position variance
%                state(n).metr(:,r) = time(n) prm.rcvs{r} meteo. parameters
%                state(n).bcpr(s,r) = time(n) prm.sats{s}-prm.rcvs{r} phase bias
%                state(n).ecls(s,1) = time(n) prm.sats{s} eclipse flags
%                state(n).erp       = time(n) earth rotation parameters
%                state(n).psun      = time(n) sum position (eci)
%                state(n).pmoon     = time(n) moon position (eci)
%                state(n).clko      = time(n) satellite clock offset to gps time
%            abort   = abort status (1:aborted,0:completed,-1:error)
%            msg     = error message
%            stats   = satellite eclipse periods
% [note]   :
% [version]: $Revision: 5 $ $Date: 06/07/08 1:16 $
%            Copyright(c) 2004-2006 by T.Takasu, all rights reserved
% [history]: 06/03/26  0.1  separated from gpsestd.m
%            06/06/24  0.2  add argout msg
%-------------------------------------------------------------------------------
if nargin<5, opt=0; end, abort=0; msg=''; stats=[];
if ~abort&opt~=2, [state,abort,msg,stats]=reads(td,time,state,prm); end
if ~abort&opt~=1, [state,abort,msg]=readr(td,time,state,prm); end

% read satellite states --------------------------------------------------------
function [state,abort,msg,stats]=reads(td,time,state,prm)
abort=0; stats=[]; msg='';
srceph={'','igs','igr','igu','ephf','ephb','ephfb','ephs','igp','cod','emr','esa','gfz','jpl','mit','codr'};
srcclk={'','igs','igr','igu','clkf','clkb','clkfb','igp','cod','emr','esa','gfz','jpl','mit','igscod','codr','igrcodr'};
srcerp={'','igs','igr','bulb','bula','c04','erpf','erpb','erpfb','cod','emr','esa','gfz','jpl','mit'};
nt=length(time); ns=length(prm.sats);

s.ephs=[]; s.ephv=[]; s.clks=[]; s.clkv=[];
if any(prm.est.sato>=2), s.ephs=repmat(nan,6,ns); end
if any(prm.est.satc>=2), s.clks=repmat(nan,1,ns); end
if prm.sdfact2>0
    if any(prm.est.sato>=2), s.ephv=repmat(nan,6,ns); end
    if any(prm.est.satc>=2), s.clkv=repmat(nan,1,ns); end
end
s.erp=prm.erp; s.ecls=uint8(zeros(ns,1)); s.clko=0;
state=repmat(s,nt,1);

% satellite ephemerides/clock from navigation message
if gmsg('reading navigation msgs : %s %s',pstr(td,time),prm.src.nav)
    abort=1; return;
end
[nav,inav]=readnav(td,time,prm.sats,prm.rcvs,prm.dirs.nav,prm.src.nav);
if isempty(nav)
    gt_log('no navigation message   : %s dir=%s',prm.src.nav,prm.dirs.nav);
    abort=-1; msg='no navigation message'; return;
end
for n=1:ns
    navs=nav(inav==n,:);
    for m=1:nt
        [pos,dts,vel]=navtostate(td,time(m),navs);
        state(m).ephs(:,n)=[pos;vel];
        state(m).clks(:,n)=dts(1);
    end
end
for n=2:length(srceph) % satellite ephemerides
    i=find(prm.est.sato==n);
    if ~isempty(i)
        if gmsg('reading satellite ephems : %s %s',pstr(td,time),srceph{n})
            abort=1; return;
        end
        if 5<=n&n<=7, dirs=prm.dirs.inp; else dirs=prm.dirs.eph; end
        [eph,sig]=readeph(td,time,prm.sats(i),dirs,srceph{n},prm.tunit,'interp');
        for m=1:nt
            state(m).ephs(:,i)=eph(m,:,:);
            if prm.sdfact2>0, state(m).ephv(:,i)=sig(m,:,:).^2; end
        end
    end
end
for n=2:length(srcclk) % satellite clock
    i=find(prm.est.satc==n);
    if ~isempty(i)
        if gmsg('reading satellite clocks : %s %s',pstr(td,time),srcclk{n})
            abort=1; return;
        end
        if 5<=n&n<=7, dirs=prm.dirs.inp; else dirs=prm.dirs.clk; end
        if prm.clkintp
            [clk,sig]=readclk(td,time,prm.sats(i),dirs,srcclk{n},prm.tuinp,'interp');
        else
            [clk,sig]=readclk(td,time,prm.sats(i),dirs,srcclk{n},prm.tuinp);
        end
        for m=1:nt
            c=squeeze(clk(m,1,:))';
            j=find(~isnan(c)&~isnan(state(m).clks(1,:)));
            if ~isempty(j), state(m).clko=mean(c(j)-state(m).clks(1,j)); end
            state(m).clks(:,i)=clk(m,:,:);
            if prm.sdfact2>0, state(m).clkv(:,i)=sig(m,:,:).^2; end
        end
    end
end
if prm.est.erp>=2 % earth rotation parameters
    if gmsg('reading erp : %s %s',pstr(td,time),srcerp{prm.est.erp})
        abort=1; return;
    end
    for m=1:nt
        utc_tai=prm_utc_tai(td+time(m)/86400,1);
        tu=td+(time(m)+19+utc_tai)/86400;
        state(m).erp=readerp(tu,prm.dirs.erp,srcerp{prm.est.erp});
        if prm.erpvar
            state(m).erp(1:3)=state(m).erp(1:3)+erpvar(tu,utc_tai);
        end
    end
end
if gmsg('computing satellite eclipse : %s',pstr(td,time)), abort=1; return; end
t=repmat(-inf,ns,1);
for m=1:nt
    utc_tai=prm_utc_tai(td+time(m)/86400,1);
    tu=td+(time(m)+19+utc_tai)/86400;
    
    % eci to ecef transformation matrix
    [U,q,q,q,q,q,du]=ecsftoecef(tu,state(m).erp,utc_tai,prm.nutmodel);
    
    % sun/moon positions
    [state(m).rsun,state(m).rmoon]=sunmoonpos(tu);
    
    % satellite eclipse conditions (eci)
    psun=U*state(m).rsun; pmoon=U*state(m).rmoon;
    for n=1:ns
        pos=state(m).ephs(1:3,n);
        if all(~isnan(pos))&shadowfunc(pos,psun,pmoon)<1
            state(m).ecls(n)=1; t(n)=time(m);
        elseif time(m)<=t(n)+prm.ecltime % post-eclipse manuever
            state(m).ecls(n)=2;
        end
    end
    % transform satellite ecef position/velocity to eci
    if prm.est.erp~=1
        [state(m).ephs,state(m).ephv]=toeci(state(m).ephs,state(m).ephv,U,du);
    end
end
ecls=[state.ecls];
for n=1:ns
    i=find(ecls(n,1:end-1)~=1&ecls(n,2:end)==1)+1;
    j=find(ecls(n,1:end-1)==1&ecls(n,2:end)~=1);
    if ecls(n,1)==1, i=[1,i]; end, if ecls(n,end)==1, j=[j,nt]; end
    stats=[stats;repmat(n,length(i),1),time(i)',time(j)'];
end

% transform satellite ecef position/velocity to eci ----------------------------
function [ephs,ephv]=toeci(ephs,ephv,U,du)
ephs(4:6,:)=U'*ephs(4:6,:)+du'*ephs(1:3,:);
ephs(1:3,:)=U'*ephs(1:3,:);
if ~isempty(ephv)
    ephv(1:3,:)=(U'*sqrt(ephv(1:3,:))).^2;
    ephv(4:6,:)=(U'*sqrt(ephv(4:6,:))).^2;
end

% read station states ----------------------------------------------------------
function [state,abort,msg]=readr(td,time,state,prm)
abort=0; msg='';
srcclk={'','igs','igr','igu','clkf','clkb','clkfb','igp','cod','emr','esa','gfz','jpl','mit','igscod'};
srcpos={'','igssnx','itrf2000','itrf97','gsipos','posf','posb','posfb','igs00','igb00'};
srczpd={'','igssnx','igsmon','zpdf','zpdb','zpdfb','mso'};
nt=length(time); ns=length(prm.sats); nr=length(prm.rcvs);

for n=1:nt
    state(n).clkr=repmat(nan,1,nr);
    state(n).posr=repmat(nan,6,nr);
    if any(prm.est.rcvz>=2), state(n).zpdr=repmat(nan,1,nr); else state(n).zpdr=[]; end
    if any(prm.est.rcvp>=2), state(n).posv=repmat(nan,3,nr); else state(n).posv=[]; end
    if any(prm.est.rcvb>=1), state(n).bcpr=repmat(nan,nr,ns); else state(n).bcpr=[]; end
    if ~isempty(prm.metsrc), state(n).metr=repmat(nan,3,nr); else state(n).metr=[]; end
end
for n=2:length(srcclk) % receiver clock
    i=find(prm.est.rcvc==n);
    if ~isempty(i)
        if gmsg('reading receiver clocks : %s %s',pstr(td,time),srcclk{n})
            abort=1; return
        end
        if 5<=n&n<=7, dirs=prm.dirs.inpr; else dirs=prm.dirs.clk; end
        if prm.clkintp
            [clk,sig]=readclk(td,time,prm.rcvs(i),dirs,srcclk{n},prm.tunit,'interp');
        else
            [clk,sig]=readclk(td,time,prm.rcvs(i),dirs,srcclk{n},prm.tunit);
        end
        for m=1:nt
            state(m).clkr(:,i)=clk(m,:,:);
        end
    end
end
for n=2:length(srczpd) % station tropospheric zenith path delay
    i=find(prm.est.rcvz==n);
    if ~isempty(i)
        if gmsg('reading tropos parameters : %s %s',pstr(td,time),srczpd{n})
            abort=1; return
        end
        if 4<=n&n<=6, dirs=prm.dirs.inpr; else dirs=prm.dirs.trop; end
        [zpd,sig]=readtrop(td,time,prm.rcvs(i),dirs,srczpd{n},prm.tunit,'interp');
        for m=1:nt
            state(m).zpdr(:,i)=zpd(m,1,:);
        end
    end
end
for n=2:length(srcpos) % station position
    i=find(prm.est.rcvp==n);
    if ~isempty(i)
        if gmsg('reading station positions : %s %s',pstr(td,time),srcpos{n})
            abort=1; return
        end
        if 6<=n&n<=8, dirs=prm.dirs.inpr; else dirs=prm.dirs.pos; end
        if isempty(prm.tp)
            [pos,cov]=readpos(td,time,prm.rcvs(i),dirs,srcpos{n},prm.tunit,'interp');
            for m=1:nt
                state(m).posr(1:3,i)=pos(m,:,:);
                state(m).posv(1:3,i)=cov(m,:,:);
            end
        else
            [pos,cov]=readpos(prm.tp,0,prm.rcvs(i),dirs,srcpos{n},prm.tunit,'interp');
            for m=1:nt
                state(m).posr(1:3,i)=pos(1,:,:);
                state(m).posv(1:3,i)=cov(1,:,:);
            end
        end
    end
end
for n=1:length(prm.rcvs) % phase bias
    ext={'f','b','fb'};
    if prm.est.rcvb(n)>=1
        if gmsg('reading phase biases : %s bcp%s',pstr(td,time),ext{prm.est.rcvb(n)})
            abort=1; return
        end
        [t,xs,vs,p,arc]=readest(td,time,'bcp',prm.rcvs{n},prm.dirs.inpr,...
                                ext{prm.est.rcvb(n)},prm.tunit);
        bcpr=repmat(nan,ns,nt);
        for m=1:size(arc,1)
            i=find(strcmp(prm.sats,p.sats{arc(m,3)}));
            if ~isempty(i)
                if arc(m,2)==t(end), arc(m,2)=arc(m,2)+p.tint; end
                j=find(arc(m,1)<=t&t<=arc(m,2));
                k=find(arc(m,1)<=time&time<=arc(m,2));
                if length(j)>1&~isempty(k)
                    bcpr(i,k)=interp1(t(j),xs(j,arc(m,3)),time(k),'linear','extrap');
                end
            end
        end
        for m=1:nt, state(m).bcpr(n,:)=bcpr(:,m)'; end
    end
end
for n=1:length(prm.rcvs) % input position from parameter file
    if prm.est.rcvp(n)>=length(srcpos)&~isempty(prm.rcvposs)
        i=find(strcmp(prm.rcvs{n},prm.rcvposs(:,1)));
        if ~isempty(i)
            pos=[prm.rcvposs{i,2:7}]';
            for m=1:nt
                state(m).posr(1:3,n)=pos(1:3,:);
                state(m).posv(1:3,n)=pos(4:6,:).^2;
            end
        else
            gt_log('no receiver position    : %s file=%s',prm.rcvs{n},prm.rcvposf);
        end
    end
end
if ~isempty(prm.metsrc) % station meteorological parameters
    if gmsg('reading meteo parameters : %s %s',pstr(td,time),prm.metsrc)
        abort=1; return
    end
    metr=readmet(td,time,prm.rcvs,prm.dirs.trop,prm.metsrc);
    for m=1:nt, state(m).metr=metr(:,:,m)'; end
end

% time string ------------------------------------------------------------------
function s=tstr(td,t), s=sprintf('%04d/%02d/%02d %02d:%02d:%02.0f',mjdtocal(td,t));
function s=pstr(td,time), s=sprintf('%s-%s',tstr(td,time(1)),tstr(td,time(end)));
