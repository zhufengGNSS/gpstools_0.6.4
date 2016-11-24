function GenRef(varargin)
%-------------------------------------------------------------------------------
% [system] : 測位システム設計解析ツール
% [module] : 参照データ生成
% [func]   : IGS PRODUCTSを読み込み参照データを生成する。
% [argin]  : (opts) = オプション
%                'td',td     = 日付(MJD)
%                'time',time = 時刻(sec)(GPST)
%                'sats',sats = 衛星
%                'rcvs',rcvs = 観測局
%                'dirs',dirs = ディレクトリ
%                'utc_tai',utc_tai = UTC-TAI(sec)
%                'tunit',tunit = 出力時間単位(sec)
% [argout] : 出力ファイル
%            ref_YYYYMMDDHH.mat : 参照データ
% [note]   : 出力は900sec単位
% [version]: $Revision: 7 $ $Date: 04/07/11 11:04 $
% [history]: 04/06/03  0.1  新規作成
%-------------------------------------------------------------------------------
p=prm_gpsest;
td=p.td; time=p.time; sats=p.sats; rcvs=p.rcvs; dirs=p.dirs; utc_tai=p.utc_tai;
tunit=p.tunit; n=1;
while n<=length(varargin)
    switch varargin{n}
    case 'td', td=varargin{n+1}; n=n+2;
    case 'time', time=varargin{n+1}; n=n+2;
    case 'sats', sats=varargin{n+1}; if ischar(sats), sats={sats}; end, n=n+2;
    case 'rcvs', rcvs=varargin{n+1}; if ischar(rcvs), sats={rcvs}; end, n=n+2;
    case 'dirs', dirs=varargin{n+1}; n=n+2;
    case 'utc_tai', utc_tai=varargin{n+1}; n=n+2;
    case 'tunit', tunit=varargin{n+1}; n=n+2;
    end
end

ts=floor(time(1)/tunit)*tunit;
te=floor((time(end)-1)/tunit)*tunit;

for time=ts:tunit:te
    GenRefData(td,time:900:time+tunit-900,sats,rcvs,dirs,utc_tai,p.odisp,p.ophas);
end

% 参照データ生成 ---------------------------------------------------------------
function GenRefData(td,time,sats,rcvs,dirs,utc_tai,odisp,ophas)

for n=1:length(time)
    tutc=td+(time(n)+19+utc_tai)/86400;
    erp(n,:)=ReadErp(tutc,dirs.erp)+[ErpVar(tutc,utc_tai),0,0];
    [U(:,:,n),P,N,gmst(n)]=EcsfToEcef(tutc,erp(n,:),utc_tai);
end
for n=1:length(sats)
    eph(:,:,n)=ReadEph(td,time,sats{n},U,dirs.eph);
    dts(:,n)=ReadClk(td,time,sats{n},dirs.clk);
end
posr=shiftdim(ReadPos(td,time(1),rcvs,dirs.pos))';
for n=1:length(rcvs)
    gpos(n,:)=EcefToGeod(posr(n,:)');
    [ap,ec]=ReadRcv(td,time(1),rcvs{n},dirs.pos);
    apc1(n,:)=ap(:,1)';
    apc2(n,:)=ap(:,2)';
    ecc(n,:)=ec';
    dtr(:,n)=ReadClk(td,time,rcvs{n},dirs.clk);
    zpd(:,n)=ReadTrop(td,time,rcvs{n},dirs.trop);
    if isnan(zpd(1,n)), zpd(:,n)=trop_saast(td,[0,pi/2],gpos(n,:)); end
end
pr=zeros(length(time),2,length(sats),length(rcvs)); azel=pr;
range=zeros(length(time),length(sats),length(rcvs)); ion=range; trop=range;
apcs=range; apcr=pr; rels=range; phw=range;
for n=1:length(rcvs)
    for m=1:length(sats)
        [pr(:,:,m,n),azel(:,:,m,n),range(:,m,n),ion(:,m,n),trop(:,m,n),...
         apcs(:,m,n),apcr(:,:,m,n),rels(:,m,n),phw(:,m,n)]=...
            RangeModel(td,time,eph(:,:,m),dts(:,m),dtr(:,n),zpd(:,n),U,gmst,...
                erp,posr(n,:)',gpos(n,:),apc1(n,:)',apc2(n,:)',ecc(n,:)',...
                utc_tai,sats,rcvs,m,n,dirs,odisp(:,:,n),ophas(:,:,n));
    end
end
epoch=MjdToCal(td,time(1));
[td,ts]=CalToMjd(epoch);
time=time(:)+(ts-time(1));
file=fullfile(dirs.ref,sprintf('ref_%04d%02d%02d%02d.mat',epoch(1:4)));
disp(['save : ',file])
save(file,'td','time','erp','eph','dts','dtr','zpd','pr','azel','range','ion',...
     'trop','apcs','apcr','rels','phw','sats','rcvs')

% 測距モデル -------------------------------------------------------------------
function [pr,azel,range,ion,trop,apcs,apcr,rels,phw]=...
    RangeModel(td,time,eph,dts,dtr,zpd,U,gmst,erp,posr,gpos,apc1,apc2,ecc,...
               utc_tai,sats,rcvs,isat,ircv,dirs,odisp,ophas)
persistent phs
if size(phs,1)~=length(sats)|size(phs,2)~=length(rcvs)
    phs=zeros(length(sats),length(rcvs));
end
C=299792458; f1=1.57542E9; f2=1.22760E9;
pr=repmat(nan,length(time),2); azel=repmat(nan,length(time),2);
range=repmat(nan,length(time),1); ion=range; trop=range; apcs=range; apcr=azel;
rels=range; phw=range;
for n=1:length(time)
    tutc=td+(time(n)+19+utc_tai)/86400;
    [rsun,rmoon]=SunMoonPos(tutc);
    posd=posr+SiteDisp(tutc,posr,U(:,:,n)*rsun,U(:,:,n)*rmoon,odisp,ophas,...
                       gmst(n),erp(n,:));
    [rsat,range(n)]=SatRange(eph(n,:)',posd,U(:,:,n),dtr(n));
    azel(n,:)=SatAzEl(U(:,:,n)*rsat,posd);
    if azel(n,2)>0
        rrcv=U(:,:,n)'*posd;
        ion(n)=ion_tec(td,time(n),azel(n,:),gpos,dirs.ion);
        apcs(n)=SatApc(rsat,rrcv,rsun,sats{isat});
        apcr(n,:)=RcvApc(azel(n,:),apc1,apc2,ecc);
        trop(n)=mapf_nmf(td,azel(n,:),gpos)*zpd(n);
        rels(n)=RelCorr(rsat,eph(n,4:6)',rrcv);
        phs(isat,ircv)=PhWindup(rsat,rsun,posd,U(:,:,n),phs(isat,ircv));
        phw(n)=phs(isat,ircv)/(2*pi);
        pr(n,:)=[range(n)+C*(dtr(n)-dts(n))+trop(n)-apcs(n)-apcr(n,1)+rels(n)+phw(n)*C/f1,...
                 range(n)+C*(dtr(n)-dts(n))+trop(n)-apcs(n)-apcr(n,2)+rels(n)+phw(n)*C/f2];
    end
end
