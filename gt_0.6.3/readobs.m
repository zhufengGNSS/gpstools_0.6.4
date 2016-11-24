function [z,iz,rpos,adel,atype,rtype,rstat,azel,slip,arc]=readobs(td,time,sats,...
             rcvs,obsdir,obssrc,carrs,ttol,tunit,ch)
%-------------------------------------------------------------------------------
% [system] : GpsTools
% [module] : read observation data
% [func]   : read observation data
% [argin]  : td   = day (mjd)
%            time = time int day (sec)
%            sats = satellites list
%            rcvs = stations list
%           (obsdir) = obs. data directory (none/'':current)
%           (obssrc) = obs. data type (none/'':'rinex')
%                 'rinex' = rinex (24hr)
%                 'rinex3'= 3-hours rinex (3hr)
%                 'rinex1'= hourly rinex (1hr)
%                 'rinexh'= high-rate rinex (15min)
%                 'obsc'  = pre-processed obs. data
%                 'obsg'  = simulated obs. data
%           (carrs)  = carriers (none/'':{'L1','L2'})
%           (ttol)   = time tolerance (sec) (inf:all data,none/empty:0.01)
%           (tunit)  = time unit to file separation (hr)
%           (ch)     = data type index
% [argout] : z  = observation data      [cp1,cp2,...,pr1,pr2,...;...]
%            iz = obervation data index [t,isat,ircv,arcf;...]
%                 t    = observation time in day (sec)
%                 isat = satellite number
%                 ircv = station number
%                 arcf = arc flag (0=na/arc cont,1=arc start,2=arc end)
%            rpos  = approx. station position (m)
%            adel  = receiver antenna delta(up/east/north) (m) 
%            atype = receiver antenna type list
%            rtype = receiver model type list
%            rstat = receiver states [t,posx,posy,posz,velx,vely,velz;...] (ecef)
%            azel  = satellite azimath/elevation angle [az,el;...] (rad)
%            slip  = cycle-slip position [t,isat,ircv,type;...]
%            arc   = arc info [ts,te,sat,rcv,n1,n2;...]
% [note]   :
% [version]: $Revision: 32 $ $Date: 06/07/21 14:06 $
%            Copyright(c) 2004-2006 by T.Takasu, all rights reserved
% [history]: 04/11/12   0.1  new
%            04/11/16   0.2  add argout rtype
%            04/11/18   0.3  check non-cc receivers
%            04/12/08   0.4  argin tunit in hour instead of sec
%            05/03/10   0.5  support 1h,3h rinex
%            05/04/25   0.6  support compressed obs. data file
%                            support p1c1 bias time dependency
%                            use C1 instead of P1 if no P1 exists
%            05/08/12   0.7  add argout rstat,azel,slip
%            06/02/14   0.8  add argin ch
%-------------------------------------------------------------------------------
if nargin<5, obsdir=''; end
if nargin<6, obssrc=''; end
if nargin<7, carrs={}; end
if nargin<8, ttol=[]; end
if nargin<9, tunit=[]; end
if nargin<10, ch=[]; end
if isempty(obssrc), obssrc='rinex'; end
if isempty(carrs), carrs={'L1','L2'}; end
if ischar(sats), sats={sats}; end
if ischar(rcvs), rcvs={rcvs}; end
if ischar(carrs), carrs={carrs}; end
if isempty(ttol), ttol=0.01; end

z=[]; iz=[]; rpos=repmat(nan,3,length(rcvs)); adel=repmat(nan,3,length(rcvs));
atype={}; rtype={}; rstat={}; azel=[]; slip=[]; arc=[];
for n=1:length(rcvs)
    [zn,izn,rp,ad,at,rt,rs,an,sl,arcn]=...
        ReadObsFiles(td,time,obsdir,obssrc,sats,rcvs{n},carrs,ttol,tunit,ch);
    if ~isempty(zn), z=[z;zn]; izn(:,3)=n; iz=[iz;izn]; azel=[azel;an]; end
    if ~isempty(rp), rpos(:,n)=rp; end
    if ~isempty(ad), adel(:,n)=ad; end
    atype{n}=at;
    rtype{n}=rt;
    rstat{n}=rs;
    if ~isempty(sl), slip=[slip;sl(:,1:2),repmat(n,size(sl,1),1),sl(:,3)]; end
    if ~isempty(arcn), arcn(:,4)=n; arc=[arc;arcn]; end
end
if ~isempty(iz), [iz,i]=sortrows(iz,[1,3,2]); z=z(i,:); end

% read obs. files --------------------------------------------------------------
function [z,iz,rpos,adel,atype,rtype,rstat,azel,slip,arc]=ReadObsFiles(td,time,...
             obsdir,obssrc,osats,rcv,carrs,ttol,tunit,ch)
z=[]; iz=[]; n=1; te=time(1); rpos=[]; adel=[]; atype=''; rtype=''; rstat=[];
azel=[]; slip=[]; arc=[];
while ~isempty(n)
    [zn,izn,te,rposn,adeln,atypen,rtypen,rstatn,azeln,slipn,arcn]=...
        ReadObsFile(td,time(n),obsdir,obssrc,osats,rcv,carrs,tunit);
    if ~isempty(zn)
        rpos=rposn; adel=adeln; atype=atypen; rtype=rtypen;
        if ~isempty(ch), zn=zn(:,ch); end
    end
    z=[z;zn]; iz=[iz;izn]; rstat=[rstat;rstatn]; slip=[slip;slipn]; arc=[arc;arcn];
    if ~isempty(azeln), azel=[azel;azeln];
    else azel=[azel;repmat(nan,size(izn,1),2)]; end
    n=min(find(te<=time));
end
if ~isempty(z)
    [iz,i]=unique(iz,'rows'); z=z(i,:); % delete duplicated data
    
    if ~isinf(ttol) % extract specified time
        index=zeros(size(iz,1),1); i=1;
        for n=1:length(time)
            for i=i:size(iz,1)
                if time(n)-ttol<=iz(i,1)&iz(i,1)<=time(n)+ttol, index(i)=1;
                elseif time(n)<iz(i,1), break, end
            end
        end
        i=find(index); z=z(i,:); iz=iz(i,:);
    end
end

% read obs. file ---------------------------------------------------------------
function [z,iz,te,rpos,adel,atype,rtype,rstat,azel,slip,arc]=...
             ReadObsFile(td,ts,obsdir,obssrc,osats,rcv,carrs,tunit)
z=[]; iz=[]; te=ts+86400; epoch=[]; data=[]; index=[]; time=[]; rpos=[]; adel=[];
atype=''; rtype=''; rstat=[]; azel=[]; slip=[]; arc=[];
dt=mjdtocal(td,ts); tdd=caltomjd(dt(1:3));
if length(rcv)>=4, rcvf=rcv(end-3:end); else rcvf=rcv; end
switch obssrc,
case 'rinex' % rinex (24hr)
    day=tdd-caltomjd([dt(1),1,1])+1;
    file=sprintf('%s%03d%1d.%02do',rcvf,day,0,mod(dt(1),100));
    file=gfilepath(obsdir,file,dt,rcvf);
    [epoch,time,types,units,sats,rcv,data,index,rpos,adel,atype,rtype]=...
        readrinexobs_(file);
    if ~isempty(epoch)
        te=(tdd-td)*86400+86400; %gmsg([],['load : ',file]);
        index=[index,zeros(size(index,1),2)];
    end

case 'rinex3' % 3-hour rinex (3hr)
    day=tdd-caltomjd([dt(1),1,1])+1;
    n=floor(dt(4)/3);
    file=sprintf('%s%03d%1d.%02do',rcvf,day,n+1,mod(dt(1),100));
    file=gfilepath(obsdir,file,dt,rcvf);
    [epoch,time,types,units,sats,rcv,data,index,rpos,adel,atype,rtype]=...
        readrinexobs_(file);
    if ~isempty(epoch)
        te=(tdd-td)*86400+10800*(n+1); %gmsg([],['load : ',fils]);
        index=[index,zeros(size(index,1),2)];
    end

case 'rinex1' % hourly rinex (1hr)
    day=tdd-caltomjd([dt(1),1,1])+1;
    file=sprintf('%s%03d%c.%02do',rcvf,day,'a'+dt(4),mod(dt(1),100));
    file=gfilepath(obsdir,file,dt,rcvf);
    [epoch,time,types,units,sats,rcv,data,index,rpos,adel,atype,rtype]=...
        readrinexobs_(file);
    if ~isempty(epoch)
        te=(tdd-td)*86400+3600*dt(5); %gmsg([],['load : ',fils]);
        index=[index,zeros(size(index,1),2)];
    end

case 'rinexh' % high-rate rinex (15min)
    day=tdd-caltomjd([dt(1),1,1])+1;
    n=floor(dt(5)/15);
    file=sprintf('%s%03d%c%02d.%02do',rcvf,day,'a'+dt(4),n*15,mod(dt(1),100));
    file=gfilepath(obsdir,file,dt,rcvf);
    [epoch,time,types,units,sats,rcv,data,index,rpos,adel,atype,rtype]=...
        readrinexobs_(file);
    if ~isempty(epoch)
        te=(tdd-td)*86400+dt(4)*3600+900*(n+1); %gmsg([],['load : ',file]);
        index=[index,zeros(size(index,1),2)];
    end

case {'obsc','obsg'}, % pre-processed/simulated obs. data
    if isempty(tunit), tunit=24; end
    n=floor(dt(4)/tunit);
    file=sprintf('%s_%s_%04d%02d%02d%02d.mat',obssrc,rcv,dt(1:3),tunit*n);
    file=gfilepath(obsdir,file,dt,rcvf);
    if exist(file), load(file), end %gmsg([],['load : ',file]);
    te=(tdd-td)*86400+tunit*3600*(n+1);
end
if ~isempty(epoch)&~isempty(data)&~isempty(index)&~isempty(time)
    [tdd,tss]=caltomjd(epoch);
    toff=(tdd-td)*86400+tss;
    index=SatIndex(index,sats,osats); i=find(index(:,1)>0);
    if ~isempty(i)
        iz=[time(i)+toff,index(i,:)];
        if any(strcmp(obssrc,{'rinex','rinex3','rinex1','rinexh'}))
            z=ArrangeObsData(data(i,:),types,carrs,iz,osats,rtype,tdd);
        else
            z=data(i,:);
        end
    end
    if ~isempty(rstat), rstat(:,1)=rstat(:,1)+toff; end
    if ~isempty(slip)
        slip(:,2)=SatIndex(slip(:,2),sats,osats); i=find(slip(:,2)>0);
        slip=[slip(i,1)+toff,slip(i,2:end)];
    end
    if ~isempty(arc)
        arc(:,3)=SatIndex(arc(:,3),sats,osats); i=find(arc(:,3)>0);
        arc=[arc(i,1:2)+toff,arc(i,3:end)];
    end
end

% read rinex obs file ----------------------------------------------------------
function [epoch,time,types,units,sats,rcv,data,index,rpos,adel,atype,rtype]=...
            readrinexobs_(file)
epoch=[]; time=[]; types={}; units={}; sats={}; rcv=''; data=[]; index=[];
rpos=[]; adel=[]; atype=''; rtype='';

cf=file; cf(end)='d'; org='';
fs={file,cf,[file,'.Z'],[file,'.gz'],[cf,'.Z'],[cf,'.gz']};
for n=1:length(fs), if exist(fs{n}), org=fs{n}; break, end, end
if ~isempty(org)
    [file,stat]=uncompact(org,'org','dir',pwd);
    if stat==0
        [epoch,time,types,units,sats,rcv,data,index,rpos,adel,atype,rtype]=...
            readrinexobs(file);
        if ~strcmp(file,org)&exist(file), delete(file); end
    end
end
if ~isempty(data)
    i=find(all(isnan(data))); if ~isempty(i), types(i)=[]; data(:,i)=[]; end
end

% convert satellite indexes ----------------------------------------------------
function indexx=SatIndex(index,sats,osats)
indexx=index; indexx(:,1)=0;
if isempty(indexx), return, end
for n=1:length(osats)
    isat=min(find(strcmp(sats,osats{n})));
    if ~isempty(isat), indexx(find(index(:,1)==isat),1)=n; end
end

% arrange obs. data by types ---------------------------------------------------
function z=ArrangeObsData(data,types,carrs,iz,sats,rtype,td)
nc=length(carrs);
z=repmat(nan,size(data,1),nc*2);
for n=1:nc % phase
    i=min(find(strcmp(types,carrs{n})));
    if ~isempty(i), z(:,n)=data(:,i); end
end
jc1=min(find(strcmp(types,'C1')));
jp1=min(find(strcmp(types,'P1')));
jp2=min(find(strcmp(types,'P2')));
for n=1:nc % code
    switch carrs{n}
    case 'L1'
        if ~isempty(jp1)
            z(:,nc+n)=data(:,jp1); % p1
        elseif ~isempty(jc1),
            z(:,nc+n)=CorrBiasP1C1(data(:,jc1),iz,sats,td); % c1+f(i)
        end
    case 'L2'
        if ~isempty(jp1)&~isempty(jp2)
            z(:,nc+n)=data(:,jp2); % p2
        elseif ~isempty(jc1)&~isempty(jp2)
            if IsNonCC(rtype)
                z(:,nc+n)=data(:,jp2); % p2 (non-cc)
            else
                z(:,nc+n)=CorrBiasP1C1(data(:,jp2),iz,sats,td); % p2+f(i) (cc)
            end
        end
    otherwise
        i=min(find(strcmp(types,['P',carrs{n}(2)])));
        if ~isempty(i), z(:,nc+n)=data(:,j); end
    end
end

% correct p1-c1 bias -----------------------------------------------------------
function z=CorrBiasP1C1(z,iz,sats,td)
C=299792458;
[dirs,f]=fileparts(which(mfilename));
dirs=fullfile(dirs,'data'); file=fullfile(dirs,'dcbs_p1c1.m');
if exist(file)
    wd=cd; cd(dirs); bias=feval('dcbs_p1c1'); cd(wd);
    bias(bias==-9.999E9)=0;
else
    gt_log('no p1-c1 dcb parameter  : %s',file); return;
end
i=max(find(bias(:,1)<=td)); if isempty(i), return, end
for n=1:length(sats)
    j=find(iz(:,2)==n); k=sscanf(sats{n},'GPS%d');
    if ~isempty(j)&~isempty(k)
        z(j,:)=z(j,:)+C*bias(i,k+1)*1E-9;
    end
end

% check non-cc receivers reporting C1 instead of P1 ----------------------------
function f=IsNonCC(rtype)
[dirs,f]=fileparts(which(mfilename));
file=fullfile(dirs,'data','rcvs_noncc.txt');
if exist(file)
    noncc=textread(file,'%s','delimiter',',','commentstyle','matlab');
else
    gt_log('no non-cc rcv parameter : %s',file); return;
    noncc={};
end
f=0;
for n=1:length(noncc)
    if ~isempty(findstr(noncc{n},rtype)), f=1; break, end
end
