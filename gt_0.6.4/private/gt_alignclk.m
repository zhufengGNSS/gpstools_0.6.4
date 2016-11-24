function xs=gt_alignclk(td,time,xs,ix,stats,prm)
%-------------------------------------------------------------------------------
% [system] : GpsTools
% [module] : quality control and alignement of estimated clock
% [func]   : quality control and alignement of estimated clock
% [argin]  : td,time = date(mjd-gpst),time vector(sec)
%            xs,ix   = estimated states/state index
%            stats   = statistics
%            prm     = processing parameter struct
% [argout] : xs      = processed estimated clocks
% [note]   :
% [version]: $Revision: 12 $ $Date: 2008-11-25 10:02:15 +0900 (火, 25 11 2008) $
%            Copyright(c) 2004-2006 by T.Takasu, all rights reserved
% [history]: 06/04/22  0.1  new
%            06/07/04  0.2  add phase-bias correction of estimated clock
%-------------------------------------------------------------------------------
C=299792458;

% phase-bias correction of estimated clock
if prm.clkcorr
    off=0;
    for n=1:length(stats)
        if strcmp(stats{n}{1},prm.clkref), off=stats{n}{2}(6); end
    end
    for n=find(prm.est.satc==1)'
        i=ix.satc{n}(1);
        for m=1:length(stats)
            if strcmp(stats{m}{1},prm.sats{n})
                xs(:,i)=xs(:,i)-stats{m}{2}(6)-off;
            end
        end
    end
    for n=find(prm.est.rcvc==1)'
        i=ix.rcvc{n}(1);
        for m=1:length(stats)
            if strcmp(stats{m}{1},prm.rcvs{n})
                xs(:,i)=xs(:,i)+stats{m}{2}(6)-off;
            end
        end
    end
end
i=find(prm.est.satc==1);
j=[ix.satc{i}]; if prm.sclkmodel==1, j=j(1:2:end); end
if isempty(j), return, end

% quality control satellite clock compared with reference clock
if ~isempty(prm.clkqc)
    xs(:,j)=qcclk(td,time,xs(:,j)/C,prm.sats(i),prm)*C;
end
% align satellite clock to reference clock
if ~isempty(prm.clkalign)
    xs(:,j)=alignclk(td,time,xs(:,j)/C,prm.sats(i),prm)*C;
end

% phase-bias correction of estimated clock -------------------------------------
function clk=corrcbias(clk,names,stats)

% quality control clock compared with reference clock --------------------------
function clks=qcclk(td,time,clks,sats,prm)
maxdclk=0.1; % qc threshold (ns)
if strcmp(prm.clkqc,'igscod'), tint=30; else tint=300; end
t=(floor(time(1)/tint):ceil(time(end)/tint))*tint;
clkr=squeeze(readclk(td,t,sats,prm.dirs.clk,prm.clkqc));
if all(all(isnan(clkr)))
    gt_log('satellite clock missing : %s %s',pstr(td,t(1),t(end)),prm.clkqc);
    return;
end
[tt,i,j]=intersect(time,t);
for n=1:length(tt)
    k=find(~isnan(clks(i(n),:))&~isnan(clkr(j(n),:)));
    if ~isempty(k)
        dc=(clks(i(n),:)-mean(clks(i(n),k)))-(clkr(j(n),:)-mean(clkr(j(n),k)));
        for k=find(abs(dc)*1E9>maxdclk) % exclude clock
            m=find(tt(n)-tint<time&time<tt(n)+tint);
            clks(m,k)=nan;
            gt_log('satellite clock excluded : %s %s dclk=%5.2fns',sats{k},...
                   pstr(td,time(m(1)),time(m(end))),dc(k)*1E9);
        end
    end
end

% align clock to reference clock -----------------------------------------------
function clks=alignclk(td,time,clks,sats,prm)
if strcmp(prm.clkalign,'code'), tint=30; else tint=300; end
t=(floor(time(1)/tint):ceil(time(end)/tint))*tint;
clkr=squeeze(readclk(td,t,sats,prm.dirs.clk,prm.clkalign));
if all(all(isnan(clkr)))
    gt_log('satellite clock missing : %s %s',pstr(td,t(1),t(end)),prm.clkalign);
    return;
end
[tt,i,j]=intersect(time,t);
for n=1:size(clks,2)
    off=clks(:,n)-interp1(tt,clks(i,n),time(:));
    off2=clks(i(end)+1:end,n)-clks(i(end),n);
    clks(:,n)=interp1(t,clkr(:,n),time(:))+off;
    clks(i(end)+1:end,n)=clks(i(end),n)+off2;
end

% time string ------------------------------------------------------------------
function s=tstr(td,t)
s=sprintf('%04d/%02d/%02d %02d:%02d:%02.0f',mjdtocal(td,t));

function s=pstr(td,ts,te)
t=mjdtocal(td,te); s=sprintf('%s-%02d:%02d:%02.0f',tstr(td,ts),t(4:6));
