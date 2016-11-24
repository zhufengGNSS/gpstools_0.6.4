function ChkObs(varargin)
%-------------------------------------------------------------------------------
% [system] : 測位システム設計解析ツール
% [module] : 観測データチェック
% [func]   : GpsObsで生成した電離層/バイアスフリー搬送波位相の品質をチェックする。
% [argin]  : (opts) = オプション
%                'sats',sats : 衛星名
%                'rcvs',rcvs : 観測局名
%                'obs'  : 観測データ有無表示
%                'mpc'  : マルチパス補正
%                'cent' : 中央表示
%                'range',range : レンジ
% [argout] : なし
% [note]   :
% [version]: $Revision: 8 $ $Date: 04/11/09 21:06 $
% [history]: 04/06/03  0.1  新規作成
%-------------------------------------------------------------------------------
global p_, p_=prm_gpsest;
sats=p_.sats; rcvs=p_.rcvs; obs=0; mpc=0; cent=0; range=[]; n=1; 
while n<=length(varargin)
    switch varargin{n}
    case 'sats', sats=varargin{n+1}; n=n+2;
    case 'rcvs', rcvs=varargin{n+1}; n=n+2;
    case 'obs', obs=1; n=n+1;
    case 'mpc', mpc=1; n=n+1;
    case 'cent', cent=1; n=n+1;
    case 'range', range=varargin{n+1}; n=n+2;
    end
end
if ischar(sats), sats={sats}; end
if ischar(rcvs), rcvs={rcvs}; end
time=p_.time(1:end-1);

if obs, obsstat(p_.td,time,sats,rcvs,p_.dirs)
else obsref(p_.td,time,sats,rcvs,mpc,cent,range,p_.dirs,p_.elmin)
end

% 観測データ有無表示 -----------------------------------------------------------
function obsstat(td,time,sats,rcvs,dirs)

str=sprintf('%02d/%02d/%02d %02d:%02d:%02.0f-',MjdToCal(td+time(1)/86400));
for n=1:length(rcvs)
    figure, hold on, grid on, label='';
    [z,iz,zsats]=LoadObs(td,time,rcvs{n},dirs); clear z
    
    for m=1:length(sats)
        obs=repmat(nan,length(time),1);
        k=find(strcmp(sats{m},zsats));
        if ~isempty(k)
            [t,i]=intersect(time,iz(find(iz(:,2)==k),1)); obs(i)=0;
        end
        plot(time/3600,obs+100*(length(sats)-m)+50,'.')
        label=[sats{m},'|',label];
    end
    set(gca,'ytick',50:100:length(sats)*100,'yticklabel',label)
    title([rcvs{n},' 観測データ有無 : ',str])
    xlabel('time(H)'), axis([time([1,end])/3600,0,100*length(sats)])
end

% 観測データ比較 ---------------------------------------------------------------
function obsref(td,time,sats,rcvs,mpc,cent,range,dirs,elmin)
C=299792458; f1=1.57542E9; f2=1.22760E9;
lam1=C/f1; lam2=C/f2; cif=[f1^2;-f2^2]/(f1^2-f2^2);
utc_tai=-32;

str=sprintf('%02d/%02d/%02d %02d:%02d:%02d',MjdToCal(td));

[tr,prr,rr,azel,eph,rsats,rrcvs]=LoadRef(td,time,dirs); % 参照データ

for n=1:length(rcvs)
    ircv=find(strcmp(rcvs{n},rrcvs));
    if isempty(ircv), disp(['no rcv data : ',rcvs{n}]), continue, end
    
    [z,iz,zsats]=LoadObs(td,time,rcvs{n},dirs);
    if mpc, [mpcc,mpcs]=LoadMpc(td,time,rcvs{n},dirs); end
    
    ert=[]; mert=[];
    for m=1:length(sats)
        isat=find(strcmp(sats{m},rsats));
        zsat=find(strcmp(sats{m},zsats));
        if isempty(isat)|isempty(zsat)
            disp(['no sat data : ',sats{m}]), continue
        end
        i=find(iz(:,2)==zsat);
        tz=iz(i,1); zz=z(i,:); izz=iz(i,:);
        [tt,ti,tj]=intersect(tz,tr);
        azelt=azel(tj,:,isat,ircv);
        epht=eph(tj,:,isat);
        
        figure, hold on, grid on
        is=find(izz(:,4)==1)';
        ie=find(izz(:,4)==2)'; if length(ie)<length(is), ie=[ie,size(izz,1)]; end
        is=is(1:length(ie));
        for i=[is;ie]
            j=find(tz(i(1))<=tt&tt<=tz(i(2))&azelt(:,2)>=elmin);
            pr=prr(tj(j),:,isat,ircv)*cif;
            if mpc
                for k=1:length(pr)
                    mpr=RcvMpc(azelt(j(k),:),mpcc,mpcs);
                    pr(k)=pr(k)+mpr;
                end
            end
            err=zz(ti(j),5)-pr;
            erf=err(~isnan(err));
            if ~isempty(erf), merr=mean(erf); mstd=std(erf);
            else, merr=0; mstd=0; end
            if ~isempty(j)
                if cent, err=err-merr; end
                plot(tt(j)/3600,err,'-')
                plot(tt(j)/3600,err,'.')
                plot(tt(j(1))/3600,err(1),'r.')
                plot(tt(j(end))/3600,err(end),'g.')
                disp(sprintf('%s-%s : t=%6d-%6d n=%4d mean=%7.4fm std=%6.4fm',...
                     sats{m},rcvs{n},tt(j(1)),tt(j(end)),length(err),merr,mstd))
            end
            ert=[ert;erf]; mert=[mert;erf-merr];
        end
        xlabel('time(H)'), ylabel('bias error(m)')
        if isempty(range), range=[-inf,inf]; end
        axis([time(1)/3600,time(end)/3600,range])
        title([sats{m},'-',rcvs{n},' 観測データ-モデル差(IF搬送波位相) : ',str])
    end
    disp(sprintf('TOTAL %s :                 n=%4d mean=%7.4fm std=%6.4fm',...
         rcvs{n},length(ert),mean(ert),std(mert)))
end

% 観測データロード -------------------------------------------------------------
function [z,iz,sats]=LoadObs(td,time,rcv,dirs)
z=[]; iz=[]; tunit=10800;
ts=floor(time(1)/tunit)*tunit:tunit:time(end);
for n=1:length(ts)
    epoch=MjdToCal(td,ts(n));
    file=fullfile(dirs.obc,sprintf('obsc_%s_%04d%02d%02d%02d.mat',rcv,epoch(1:4)));
    if exist(file)
        load(file), disp(['load : ',file])
        [tdd,tss]=CalToMjd(epoch); time=time+(tdd-td)*86400+tss;
        z=[z;data]; index(:,2)=n; iz=[iz;time,index];
    else, disp(['no file : ',file]), end
end
[iz,i]=sortrows(iz,[1,3,2]); z=z(i,:);

% 参照データロード -------------------------------------------------------------
function [t,prr,rr,azelr,ephr,sats,rcvs]=LoadRef(tdd,time,dirs)
t=[]; prr=[]; rr=[]; apcsr=[]; apcrr=[]; ephr=[]; azelr=[]; tunit=10800;
ts=floor(time(1)/tunit)*tunit:tunit:time(end);
for n=1:length(ts)
    epoch=MjdToCal(tdd,ts(n));
    file=fullfile(dirs.ref,sprintf('ref_%04d%02d%02d%02d.mat',epoch(1:4)));
    if exist(file)
        load(file), disp(['load : ',file])
        t=[t;time+(td-tdd)*86400];
        prr=cat(1,prr,pr);
        rr=cat(1,rr,range);
        apcsr=cat(1,apcsr,apcs);
        apcrr=cat(1,apcrr,apcr);
        ephr=cat(1,ephr,eph);
        azelr=cat(1,azelr,azel);
    else, disp(['no file : ',file]), end
end
rr=[rr-apcsr-squeeze(apcrr(:,1,:,:)),rr-apcsr-squeeze(apcrr(:,2,:,:))];

% マルチパスデータロード ------------------------------------------------------
function [mpcc,mpcs]=LoadMpc(td,time,rcv,dirs)
epoch=MjdToCal(td,time(1));
file=fullfile(dirs.ref,sprintf('mpps_%s_%04d%02d%02d.mat',rcv,epoch(1:3)));
if exist(file)
    load(file), disp(['load : ',file])
else, disp(['no file : ',file]), end
