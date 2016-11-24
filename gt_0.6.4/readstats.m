function [t,index,ress,outs,sats,rcvs,prm]=readstats(td,time,dirs,fb,tunit,rcvs)
%-------------------------------------------------------------------------------
% [system] : GpsTools
% [module] : read residuals
% [func]   : read residuals
% [argin]  : td    = date (mjd-gpst)
%            time  = time vector (sec)
%           (dirs) = estimation data directory  (default:current)
%           (fb)   = estimation direction ('resf':forward,'resb':backward)
%           (tunit)= processing unit time (hr) (default:24)
%           (rcvs) = stations list
% [argout] : t     = estimation time vector (sec)
%            index = satellite/station index [sat1,sat2,rcv1,rcv2;...]
%            ress  = residuals [prefit,postfit,phase-bias;...]
%            outs  = outlier flags
%            sats  = satellites list
%            rcvs  = stations list
%            prm   = processing parameters
% [note]   :
% [version]: $Revision: 12 $ $Date: 2008-11-25 10:02:15 +0900 (火, 25 11 2008) $
%            Copyright(c) 2004-2006 by T.Takasu, all rights reserved
% [history]: 04/11/17  0.1  new
%            05/07/29  0.2  separate stats files
%            06/07/03  0.3  add argout prm
%-------------------------------------------------------------------------------
if nargin<3, dirs=''; end
if nargin<4, fb='resf'; end
if nargin<5, tunit=24; end
if nargin<6, rcvs={}; end

if isempty(rcvs)
    [t,index,ress,outs,sats,rcvs,prm]=ReadStatsC(td,time,dirs,fb,tunit);
else
    t=[]; index=[]; ress=[]; outs=[]; sats={};
    for n=1:length(rcvs)
        [t_,index_,ress_,outs_,sats_,prm_]=ReadStatsR(td,time,dirs,fb,tunit,rcvs{n});
        if ~isempty(t_)
            if isempty(sats), sats=sats_; prm=prm_; end
            index_(:,3)=n;
            t=[t;t_]; index=[index;index_]; ress=[ress;ress_]; outs=[outs;outs_];
        end
    end
end

% read ress from combined stats file -------------------------------------------
function [t,indexs,ress,outs,sats,rcvs,prm]=ReadStatsC(td,time,dirs,fb,tunit)
tu=tunit*3600;
t=[]; indexs=[]; ress=[]; outs=[]; sats={}; rcvs={}; prm=[];
for ts=floor(time(1)/tu)*tu:tu:time(end)
    ep=mjdtocal(td,ts);
    file=sprintf('%s_%04d%02d%02d%02d.mat',fb,ep(1:4));
    file=gfilepath(dirs,file,ep);
    if exist(file)
        epoch=[]; time=[]; index=[]; residual=[]; outl=[];
        load(file)
        if ~isempty(epoch)
            [tdd,tss]=caltomjd(epoch);
            t=[t;time+(tdd-td)*86400+tss];
            indexs=[indexs;index];
            ress=[ress;residual];
            if ~isempty(outl)
                outs=[outs;outl];
            else
                outs=[outs;zeros(length(time),1)];
            end
            if ~isempty(prm)
                sats=prm.sats;
                rcvs=prm.rcvs;
            end
        end
    end
end
index=indexs;

% read ress from station stats file --------------------------------------------
function [t,indexs,ress,outs,sats,prm]=ReadStatsR(td,time,dirs,fb,tunit,rcv)
tu=tunit*3600;
t=[]; indexs=[]; ress=[]; outs=[]; sats={}; prm=[];
for ts=floor(time(1)/tu)*tu:tu:time(end)
    ep=mjdtocal(td,ts);
    file=sprintf('%s_%s_%04d%02d%02d%02d.mat',fb,rcv,ep(1:4));
    file=gfilepath(dirs,file,ep);
    if exist(file)
        epoch=[]; time=[]; index=[]; residual=[]; outl=[];
        load(file)
        if ~isempty(epoch)
            [tdd,tss]=caltomjd(epoch);
            t=[t;time+(tdd-td)*86400+tss];
            indexs=[indexs;index];
            ress=[ress;residual];
            if ~isempty(outl)
                outs=[outs;outl];
            else
                outs=[outs;zeros(length(time),1)];
            end
            if ~isempty(prm), sats=prm.sats; end
        end
    end
end
