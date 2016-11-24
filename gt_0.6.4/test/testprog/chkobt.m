function ChkObt(varargin)
%-------------------------------------------------------------------------------
% [system] : 測位システム設計解析ツール
% [module] : 衛星運動モデルチェック
% [func]   : 精密暦と衛星運動モデルを比較し精度を確認する。
% [argin]  : なし
% [argout] : なし
% [note]   :
% [version]: $Revision: 12 $ $Date: 2008-11-25 10:02:15 +0900 (轣ｫ, 25 11 2008) $
% [history]: 04/06/05  0.1  新規作成
%-------------------------------------------------------------------------------
global p_, p_=prm_gpsest;
td=p_.td; time=p_.time; sats=p_.sats;
oprm.g_nmax=8;
oprm.p_plgrv='grv_sunmoon';
oprm.p_solarpr='srp_code';
oprm.p_eclipse='';
oprm.p_tidal='';
oprm.p_relativ='';
oprm.p_deltav='';
n=1;
while n<=length(varargin)
    switch varargin{n}
    case 'td', td=varargin{n+1}; n=n+2;
    case 'time', time=varargin{n+1}; n=n+2;
    case 'sats', sats=varargin{n+1}; n=n+2;
    case 'gnmax', oprm.g_nmax=varargin{n+1}; n=n+2;
    case 'plgrv', oprm.p_plgrv=varargin{n+1}; n=n+2;
    case 'solarp', oprm.p_solarpr=varargin{n+1}; n=n+2;
    case 'eclipse', oprm.p_eclipse=varargin{n+1}; n=n+2;
    case 'tidal', oprm.p_tidal=varargin{n+1}; n=n+2;
    case 'relativ', oprm.p_relativ=varargin{n+1}; n=n+2;
    case 'deltav', oprm.p_deltav=varargin{n+1}; n=n+2;
    otherwise, n=n+1;
    end
end
time=time(1):900:time(end);
if ischar(sats), sats={sats}; end

ti=sprintf('%04d/%02d/%02d %02d:%02d:%02d',MjdToCal(td));
disp(sprintf('GPS ORBIT MODEL ERROR : %s : %dH',ti,(time(end)-time(1))/3600))
disp('SAT   :           POSITION(m)                         VELOCITY(m/sec)')
disp('------------------------------------------------------------------------------------')
U=zeros(3,3,length(time));
for n=1:length(time)
    tutc=td+(time(n)+19+p_.utc_tai)/86400;
    erp(n,:)=ReadErp(tutc,p_.dirs.erp)+[ErpVar(tutc,p_.utc_tai),0,0];
    U(:,:,n)=EcsfToEcef(tutc,erp(n,:),p_.utc_tai);
end
for n=1:length(sats)
    eph=ReadEph(td,time,sats{n},U,p_.dirs.eph);
    state=state_precorbit(td+(19+p_.utc_tai)/86400,time(end)-time(1),...
                          eph(1,:)',sats{n},oprm);
    E=EcsfToSatf(eph(end,:)');
    e=eph(end,:)-state';
    err(n,:)=[e(1:3)*E',e(4:6)*E'];
    disp(sprintf('%s : %6.3f(%7.3f,%7.3f,%7.3f) %10.2E(%10.2E,%10.2E,%10.2E)',...
         sats{n},norm(err(n,1:3)),err(n,1:3),norm(err(n,4:6)),err(n,4:6)))
end
rms=sqrt(mean(err.^2));
disp(sprintf('%s : %6.3f(%7.3f,%7.3f,%7.3f) %10.2E(%10.2E,%10.2E,%10.2E)',...
     'RMSE ',norm(rms(1:3)),rms(1:3),norm(rms(4:6)),rms(4:6)))
