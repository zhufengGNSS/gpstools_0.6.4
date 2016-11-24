function [nav,inav,ionprm,dutc]=readnav(td,time,sats,rcvs,navdir,navsrc)
%-------------------------------------------------------------------------------
% [system] : GpsTools
% [module] : read navigation messages
% [func]   : read broadcast navigation messages
% [argin]  : td,time = date (mjd-gpst),time vector(sec)
%           (sats)   = satellite list (default:{'GPS01',...,'GPS32'})
%           (rcvs)   = station list   (default:{})
%           (navdir) = data directory (default:current)
%           (navsrc) = navigation msg file type (default:'rinex')
%                 'brdc'  : IGS combined rinex (24hr)
%                 'rinex' : rinex (24hr)
%                 'rinex3': 3-hours rinex (3hr)
%                 'rinex1': hourly rinex (1hr)
%                 'rinexh': high-rate rinex (15min)
% [argout] : nav  = navigation messages (sorted by toe)
%                  nav(s,:) = inav(s) satellite navigation message
%                  nav(s,1:6)=Epoch : Toc (Time of Clock) (GPST)
%                           [year(2digits),month,day,hour,minute,second]
%                  nav(s,7:9)=SV Clock bias/drift/drift-rate(sec,sec/sec,sec/sec^2)
%                  nav(s,10)=IODE Issue of Data, Ephemeris
%                  nav(s,11)=Crs                      (m)
%                  nav(s,12)=Delta n                  (rad/sec)
%                  nav(s,13)=M0                       (rad)
%                  nav(s,14)=Cuc                      (rad)
%                  nav(s,15)=e Eccentricity
%                  nav(s,16)=Cus                      (rad)
%                  nav(s,17)=sqrt(A)                  (sqrt(m))
%                  nav(s,18)=Toe (Time of Ephemeris)  (sec of GPS week)
%                  nav(s,19)=Cic                      (rad)
%                  nav(s,20)=OMEGA                    (rad)
%                  nav(s,21)=Cis                      (rad)
%                  nav(s,22)=i0                       (rad)
%                  nav(s,23)=Crc                      (m)
%                  nav(s,24)=omega                    (rad)
%                  nav(s,25)=OMEGA DOT                (rad/sec)
%                  nav(s,26)=IDOT                     (rad/sec)
%                  nav(s,27)=Codes on L2 channel
%                  nav(s,28)=GPS Week # (to go with Toe)
%                  nav(s,29)=L2 P data flag
%                  nav(s,30)=SV accuracy              (m)
%                  nav(s,31)=SV health                (MSB only)
%                  nav(s,32)=TGD                      (sec)
%                  nav(s,33)=IODC Issue of Data, Clock
%                  nav(s,34)=Transmission time of message (sec of GPS week,
%                          derived e.g. from Z-count in Hand Over Word)
%                  nav(n,35:37)=spare
%            inav = satellites index
%            ionprm = ionospheric parameters [a0,a1,a2,a3;b0,b1,b2,b3]
%            dutc   = delta utc parameters [A1,A2,T,W]
% [note]   :
% [version]: $Revision: 16 $ $Date: 2008-12-12 15:49:30 +0900 (金, 12 12 2008) $
%            Copyright(c) 2004-2006 by T.Takasu, all rights reserved
% [history]: 04/11/02  0.1  new
%            05/03/10  0.2  support directory expansion
%            05/06/27  0.3  add argin navsrc
%            06/03/22  0.4  add outlier check and exclusion
%            08/06/05  0.5  support PRN32-99 (gt_0.6.3p5)
%            08/11/21  0.6  default satellite list -> {...,'GPS32'} (gt_0.6.4)
%-------------------------------------------------------------------------------
if nargin<3, sats={}; end
if nargin<4, rcvs={}; end
if nargin<5, navdir=''; end
if nargin<6, if isempty(rcvs), navsrc='brdc'; else navsrc='rinex'; end, end
if isempty(sats)
    sats={}; for n=1:32, sats={sats{:},sprintf('GPS%02d',n)}; end
end
if ischar(sats), sats={sats}; end
if ischar(rcvs), rcvs={rcvs}; end
if strcmp(navsrc,'brdc'), rcvs={''}; end
nav=[]; inav=[];

for n=1:length(rcvs)
    ts=time(1);
    while ts<=time(end)
        [navn,inavn,ionprm,dutc,tn]=readfile(td,ts,sats,rcvs{n},navdir,navsrc);
        nav=[nav;navn]; inav=[inav;inavn]; ts=tn;
    end
end
if ~isempty(nav)
    % sorted by toe
    [s,i]=sortrows(nav,[28,18]);
    nav=nav(i,:); inav=inav(i);
    
    % quality control
    [nav,inav]=qc(nav,inav,sats);
end

% read navigation message file ------------------------------------------------
function [nav,inav,ionprm,dutc,tn]=readfile(td,ts,sats,rcv,navdir,navsrc)
persistent file nsats data index nionprm ndutc, if isempty(file), file=''; end
nav=[]; inav=[]; ionprm=[]; dutc=[];
dt=mjdtocal(td,ts); tdd=caltomjd(dt(1:3));
day=tdd-caltomjd([dt(1),1,1])+1;
if length(rcv)>4, rcvf=rcv(end-3:end); else rcvf=rcv; end
switch navsrc
case {'rinex','brdc'} % rinex (24hr)
    if strcmp(navsrc,'brdc'), rcvf='brdc'; end
    f=sprintf('%s%03d%1d.%02dn',rcvf,day,0,mod(dt(1),100));
    tn=floor(ts/86400)*86400+86400;
case 'rinex3' % rinex (3hr)
    f=sprintf('%s%03d%1d.%02dn',rcvf,day,floor(dt(4)/3)+1,mod(dt(1),100));
    tn=floor(ts/10800)*10800+10800;
case 'rinex1' % rinex (1hr)
    f=sprintf('%s%03d%c.%02dn',rcvf,day,'a'+dt(5),mod(dt(1),100));
    tn=floor(ts/3600)*3600+3600;
case 'rinexh' % high-rate rinex (15min)
    f=sprintf('%s%03d%c%02d.%02dn',rcvf,day,'a'+dt(4),floor(dt(5)/15)*15,mod(dt(1),100));
    tn=floor(ts/900)*900+900;
otherwise
    warning(['not supported navigation message type : ',navsrc]);
    tn=inf; return;
end
f=gfilepath(navdir,f,dt,rcv);
if ~strcmp(f,file)|isempty(nsats)
    [nsats,nrcv,data,index,nionprm,ndutc]=readrinexnav(f); file=f;
    if isempty(nsats), return; end
    nsats={}; for i=1:99, nsats={nsats{:},sprintf('GPS%02d',i)}; end % support PRN32-99
end
inav=zeros(size(index));
for n=1:length(sats)
    i=min(find(strcmp(sats{n},nsats)));
    if ~isempty(i), inav(index==i)=n; end
end
nav=data(inav>0,:); inav=inav(inav>0); ionprm=nionprm; dutc=ndutc;

% quality control --------------------------------------------------------------
function [nav,inav]=qc(nav,inav,sats)
ex=[];
for n=1:length(sats)
    i=find(inav==n);
    if length(i)>=4
        for j=[20,22,24] % [OMG,i,omg]
            navi=nav(i,j);
            [s,k]=max(navi); navi(k)=[];
            [s,k]=min(navi); navi(k)=[];
            
            % outlier check (over 10 sigmas)
            ex=[ex;i(abs(nav(i,j)-mean(navi))>10*std(navi))];
        end
    end
end
ex=unique(ex);
for n=1:length(ex)
    e=mjdtocal((caltomjd([1980,1,6])+nav(ex(n),28)*7),nav(ex(n),18));
    gt_log('navigation msg outlier  : %s toe=%04d/%02d/%02d %02d:%02d:%02.0f',...
           sats{n},e);
end
nav(ex,:)=[]; inav(ex)=[];
