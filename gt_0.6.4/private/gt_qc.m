function [xs,xz,xa,stats]=gt_qc(td,time,xs,ps,ix,tz,iz,res,arc,prm,pass,abort)
%-------------------------------------------------------------------------------
% [system] : GpsTools
% [module] : quality control
% [func]   : quality control
% [argin]  : td,time = date(mjd-gpst),time vector(sec)
%            xs,ps   = estimated states/variances
%            tz,iz   = observation data time/index
%            res     = residuals
%            arc     = arc info
%            prm     = processing parameter struct
%            pass    = processing pass [parameter struct
%            abort   = abort status (1:abort)
% [argout] : xs      = screened estimated states
%            xz      = excluded observation data index
%            xa      = excluded arc index
%            stats   = statistics {{sat/rcv,[obs,out,resp,resf,pbr,pba,maxo]},...}
%                      sat/rcv   : satellite/receiver
%                      obs/out   : valid obs count/outlier count
%                      resp/resf : prefit/postfit residuals rms (m)
%                      pbr/pba   : phase-bias rms/average (m)
%                      maxo      : max data outage time (sec)
% [note]   : output files
% [version]: $Revision: 12 $ $Date: 2008-11-25 10:02:15 +0900 (火, 25 11 2008) $
%            Copyright(c) 2004-2006 by T.Takasu, all rights reserved
% [history]: 06/04/05  0.1  separated from gpsestd.m
%-------------------------------------------------------------------------------
xa=[]; xz=[]; stats={}; if isempty(res), return, end
i=find(iz(:,3)>0);
v=prm.vsat; vsat=[v.pout,v.rmsf,prm.maxsatout]; vsat(vsat==0)=inf;
v=prm.vrcv; vrcv=[v.pout,v.rmsf,prm.maxrcvout]; vrcv(vrcv==0)=inf;
v=prm.varc; varc=[v.pout,v.rmsf]; varc(varc==0)=inf;
for n=1:length(prm.sats)
    j=i(iz(i,1)==n);
    for m=1:length(prm.rcvs), is{n,m}=j(iz(j,2)==m); end
    s=resstats(time,tz(j),iz(j,:),res(j,:),prm);
    stats={stats{:},{prm.sats{n},s}};
    if ~abort&any([vsat(1)*s(1),vsat(2:3)]<s([2,4,7]))
        gt_log('satellite excluded : %s n=%4d out=%3d rmsf=%.4f maxo=%.0f',...
               prm.sats{n},s([1:2,4,7]));
        xz=[xz;j]; xa=[xa;find(arc(:,3)==n)];
    end
end
for n=1:length(prm.rcvs)
    j=i(iz(i,2)==n);
    s=resstats(time,tz(j),iz(j,:),res(j,:),prm);
    stats={stats{:},{prm.rcvs{n},s}};
    if ~abort&any([vrcv(1)*s(1),vrcv(2:3)]<s([2,4,7]))
        gt_log('station excluded   : %s n=%4d out=%3d rmsf=%.4f maxo=%.0f',...
               prm.rcvs{n},s([1:2,4,7]));
        xz=[xz;j]; xa=[xa;find(arc(:,4)==n)];
    end
end
for n=1:size(arc,1) % screen arc by residuals
    a=arc(n,:); j=is{a(3),a(4)}; j=j(a(1)<=tz(j)&tz(j)<=a(2));
    s=resstats(time,tz(j),iz(j,:),res(j,:),prm);
    if ~abort&any([varc(1)*s(1),varc(2)]<s([2,4]))
        gt_log('arc excluded : %s-%s %s n=%4d out=%3d rmsf=%.4f',...
               prm.sats{a(3)},prm.rcvs{a(4)},pstr(td,a(1),a(2)),s([1:2,4]));
        xz=[xz;j]; xa=[xa;n];
    end
end
% screen estimated states by variances
if pass(2)>1&~abort
    if prm.maxsatcsig>0, xs=screensatc(td,time,xs,ps,ix,prm); end
    if prm.maxrcvcsig>0, xs=screenrcvc(td,time,xs,ps,ix,prm); end
    if prm.maxrcvpsig>0, xs=screenrcvp(td,time,xs,ps,ix,prm); end
end

% statistics by residuals ------------------------------------------------------
function stats=resstats(time,tz,iz,res,prm)
rms=[0,0,0,0]; i=find(iz(:,3)==1); 
if ~isempty(i), rms=[sqrt(mean(res(i,:).^2,1)),mean(res(i,3))]; end
maxo=max(diff([time(1);sort(tz(i));time(end)+prm.tint]))-prm.tint;
stats=[size(iz,1),sum(iz(:,3)==2),rms,maxo];

% screen satellite -------------------------------------------------------------
function xs=screensatc(td,time,xs,ps,ix,prm)
C=299792458;
for n=1:length(prm.sats)
    i=ix.satc{n};
    if ~isempty(i)
        j=find(ps(:,i(1)).^2>(C*prm.maxsatcsig*1E-9)^2);
        if ~isempty(j)
            xs(j,i)=nan;
            gt_log('sat clocks deleted      : %s sig>%6.3fns n=%2d',...
                   prm.sats{n},prm.maxsatcsig,length(j));
        end
    end
end

% screen receiver clock by variences -------------------------------------------
function xs=screenrcvc(td,time,xs,ps,ix,prm)
C=299792458;
for n=1:length(prm.rcvs)
    i=ix.rcvc{n};
    if ~isempty(i)
        j=find(ps(:,i(1)).^2>(C*prm.maxrcvcsig*1E-9)^2);
        if ~isempty(j)
            xs(j,i)=nan;
            gt_log('rcv clocks deleted      : %s sig>%6.3fns n=%2d',...
                   prm.rcvs{n},prm.maxrcvcsig,length(j));
        end
    end
end

% screen receiver position by variences ----------------------------------------
function xs=screenrcvp(td,time,xs,ps,ix,prm)
for n=1:length(prm.rcvs)
    i=ix.rcvp{n};
    if ~isempty(i)
        j=find(sum(ps(:,i(1:3)).^2,2)>prm.maxrcvpsig^2);
        if ~isempty(j)
            xs(j,i)=nan;
            gt_log('rcv positions deleted   : %s sig>%6.3fm  n=%2d',...
                   prm.rcvs{n},prm.maxrcvpsig,length(j));
        end
    end
end

% time string ------------------------------------------------------------------
function s=pstr(td,ts,te)
t1=mjdtocal(td,ts); t2=mjdtocal(td,te);
s=sprintf('%04d/%02d/%02d %02d:%02d:%02.0f-%02d:%02d:%02.0f',t1,t2(4:6));
