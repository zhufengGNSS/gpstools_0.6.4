function PlotRef(varargin)
%-------------------------------------------------------------------------------
% [system] : 測位システム設計解析ツール
% [module] : 参照データ表示
% [func]   : 参照データを表示する。
% [argin]  : (opts) = オプション
%                'td',td     : 日付(MJD)
%                'time',time : 時刻ベクトル(sec)
%                'sats',sats : 衛星指定 {sat1,sat2,...}
%                'rcvs',rcvs : 観測局指定 {rcv1,rcv2,...}
%                'dirs',dirs : ディレクトリ指定
%                <type>      : データ種別
%                    'apcr' = 受信アンテナ位相中心オフセット
%                    'apcs' = 衛星アンテナ位相中心オフセット
%                    'azel' = 方位仰角
%                    'dts'  = 衛星時計
%                    'dtr'  = 観測局時計
%                    'eph'  = 衛星位置
%                    'erp'  = 地球回転パラメータ
%                    'ion'  = 電離層遅延
%                    'phw'  = phase-windup効果
%                    'pr'   = 擬似距離
%                    'range'= レンジ
%                    'rels' = 相対論効果
%                    'trop' = 対流圏遅延
%                    'zpd'  = 観測局対流圏天頂遅延
% [argout] : なし
% [note]   :
% [version]: $Revision: 12 $ $Date: 2008-11-25 10:02:15 +0900 (轣ｫ, 25 11 2008) $
% [history]: 04/06/03  0.1  新規作成
%-------------------------------------------------------------------------------
global p_, p_=prm_gpsest;
type='eph'; td=p_.td; time=p_.time; sats={}; rcvs={};
dirs=p_.dirs; n=1;
while n<=length(varargin)
    switch varargin{n}
    case 'td', td=varargin{n+1}; n=n+2;
    case 'time', time=varargin{n+1}; n=n+2;
    case 'sats', sats=varargin{n+1}; if ischar(sats), sats={sats}; end, n=n+2;
    case 'rcvs', rcvs=varargin{n+1}; if ischar(rcvs), rcvs={rcvs}; end, n=n+2;
    otherwise, type=varargin{n}; n=n+1;
    end
end
if isempty(sats), sats={''}; end
if isempty(rcvs), rcvs={''}; end
for n=1:length(rcvs)
for m=1:length(sats)
    [td,time,data,ti,yl]=LoadRef(td,time,type,sats{m},rcvs{n},dirs,p_.tunit);
    figure
    for k=1:length(yl)
        subplot(length(yl),1,k), hold on, grid on
        plot(time/3600,data(:,k),'.')
        plot(time/3600,data(:,k),'-')
        ylabel(yl{k})
        if (k==1), title(ti), end
        axis([time(1)/3600,time(end)/3600,-inf,inf])
    end
    xlabel('time (H)')
end
end

% 参照データ読み込み -----------------------------------------------------------
function [tdd,t,data,ti,yl]=LoadRef(tdd,timer,type,sat,rcv,dirs,tunit);
C=299792458; f1=1.57542E9; f2=1.22760E9; cif=[f1^2;-f2^2]/(f1^2-f2^2);
t=[]; data=[]; ti=''; yl={};
ts=floor(timer(1)/tunit)*tunit:tunit:timer(end);
for n=1:length(ts)
    dt=MjdToCal(tdd,ts(n));
    file=fullfile(dirs.ref,sprintf('ref_%04d%02d%02d%02d.mat',dt(1:4)));
    if exist(file)
        load(file)
        t=[t;time(:)+(td-tdd)*86400];
        i=find(strcmp(sats,sat));
        j=find(strcmp(rcvs,rcv));
        datan=eval(type);
        switch type
        case 'erp',
            datan(:,1:2)=datan(:,1:2)*3600*180/pi;
            data=cat(1,data,datan);
            ti='EARTH ROTATION PARAMETER'; yl={'xp(")','yp(")','UT1-UTC(sec)'};
        case 'eph',
            data=cat(1,data,datan(:,1:3,i));
            ti=[sat,' SATELITE EPHEMERIS']; yl={'xpos(m)','ypos(m)','zpos(m)'};
        case 'dts',
            data=cat(1,data,datan(:,i)*C);
            ti=[sat,' SATELLITE CLOCK BIAS']; yl={'clk bias(m)'};
        case 'dtr',
            data=cat(1,data,datan(:,j)*C);
            ti=[rcv,' RECEIVER CLOCK BIAS']; yl={'clk bias(m)'};
        case 'zpd'
            data=cat(1,data,datan(:,j));
            ti=[rcv,' TROPOSPHERIC ZENITH PATH DELAY']; yl={'zpd(m)'};
        case 'range'
            data=cat(1,data,datan(:,i,j));
            ti=[sat,'-',rcv,' RANGE']; yl={'range(m)'};
        case 'ion'
            data=cat(1,data,datan(:,i,j));
            ti=[sat,'-',rcv,' IONOSPHERIC DELAY']; yl={'ion(m)'};
        case 'trop'
            data=cat(1,data,datan(:,i,j));
            ti=[sat,'-',rcv,' TROPOSPHERIC DELAY']; yl={'trop(m)'};
        case 'apcs'
            data=cat(1,data,datan(:,i,j));
            ti=[sat,'-',rcv,' SATELLITE ANTENNA OFFSET']; yl={'apcs(m)'};
        case 'rels'
            data=cat(1,data,datan(:,i,j));
            ti=[sat,'-',rcv,' RELATIVITY EFFECT']; yl={'rels(m)'};
        case 'phw'
            data=cat(1,data,datan(:,i,j)*[C/f1,C/f2]*cif);
            ti=[sat,'-',rcv,' PHASE WIND-UP EFFECT(IF)']; yl={'phw(m)'};
        case 'pr'
            data=cat(1,data,datan(:,:,i,j));
            ti=[sat,'-',rcv,' PSEUDO RANGE']; yl={'C1(m)','P2(m)'};
        case 'azel'
            data=cat(1,data,datan(:,:,i,j)*180/pi);
            i=find(data(:,2)<0); data(i,2)=nan;
            ti=[sat,'-',rcv,' AZIMATH/ELELATION']; yl={'az(deg)','el(deg)'};
        case 'apcr'
            data=cat(1,data,datan(:,:,i,j));
            ti=[sat,'-',rcv,' ANT PHASE CENTER OFFSET']; yl={'L1(m)','L2(m)'};
        end
    else, disp(['no file : ',file]), end
end
i=find(timer(1)<=t&t<=timer(end));
t=t(i); data=data(i,:,:,:);
