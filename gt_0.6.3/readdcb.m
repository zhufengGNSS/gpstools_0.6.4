function [dcbs,sigs]=readdcb(td,time,names,dcbdir,dcbsrc)
%-------------------------------------------------------------------------------
% [system] : GpsTools
% [module] : read dcb
% [func]   : read dcb (differencial code bias) data
% [argin]  : td    = date (mjd-gpst)
%            time  = time vector (sec)
%           (names)= satellite/staion names (default:{'GPS01','GPS02',...,'GPS31'})
%           (dcbdir) = data directory (default:current)
%           (dcbsrc) = data source (defautl:'igs')
%               'igs'  = IGS ion
% [argout] : dcbs = satellite/station dcbs (m) (P1-P2)
%                dcbs(n,m) = time(n),names{n} dcb (m)
%            sigs = satellite/station dcbs std deviations (m) (P1-P2)
%                sigs(n,m) = time(n),names{n} dcb std. deviation (m)
% [note]   :
% [version]: $Revision: 6 $ $Date: 06/07/08 14:19 $
%            Copyright(c) 2004-2006 by T.Takasu, all rights reserved
% [history]: 04/11/02  0.1  new
%            05/06/07  0.2  change argin ts->time, sat,rcv->names
%                           delete argout dcbr, add argout sigs
%-------------------------------------------------------------------------------
if nargin<3|isempty(names)
    names={}; for n=1:31, names={names{:},sprintf('GPS%02d',n)}; end
end
if nargin<4, dcbdir=''; end
if nargin<5, dcbsrc='igs'; end
dcbs=repmat(nan,length(time),length(names)); sigs=dcbs;
for n=1:length(time)
    for m=1:length(names)
        [dcbs(n,m),sigs(n,m)]=readdcbdata(td,time(n),names{m},dcbdir,dcbsrc);
    end
end

% read dcb data ----------------------------------------------------------------
function [dcb,sig]=readdcbdata(td,ts,name,dcbdir,dcbsrc)
persistent file dcbss dcbrs sats rcvs, if isempty(file), file=''; end
C=299792458; dcb=nan; sig=nan; name=upper(name);

switch dcbsrc
case 'igs',
    dt=mjdtocal(td,ts); day=caltomjd(dt(1:3))-caltomjd([dt(1),1,1])+1;
    f=gfilepath(dcbdir,sprintf('igsg%03d%1d.%02di',day,0,mod(dt(1),100)),dt);
    if ~strcmp(f,file)
        file=f;
        [epoch,time,tec,rms,lats,lons,hgts,rb,dcbss,dcbrs,sats,rcvs]=...
            readionex(file);
        if isempty(epoch), return, end
    end
    i=find(strcmp(name,sats));
    if ~isempty(i), dcb=dcbss(i)*1E-9*C;
    else
        i=find(strcmp(name,rcvs));
        if ~isempty(i), dcb=dcbrs(i)*1E-9*C; end
    end

otherwise
    disp(['dcb source error : ',dcbsrc])
end
