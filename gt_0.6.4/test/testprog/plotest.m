function PlotEst(varargin)
%-------------------------------------------------------------------------------
% [system] : 測位システム設計解析ツール
% [module] : GPS軌道推定結果表示
% [func]   : GPS軌道推定結果を表示する。
% [argin]  : (opts) = オプション
%                'td',td     : 日付(MJD)
%                'time',time : 時刻ベクトル(sec)
%                'eph'       : 位置/速度推定誤差グラフ
%                'ephf'      : 位置/速度推定誤差一覧
%                'pos'       : 位置推定誤差グラフ
%                'vel'       : 速度推定誤差グラフ
%                'srp'       : 太陽輻射圧パラメータ推定結果グラフ
%                'srpf'      : 太陽輻射圧パラメータ推定結果一覧
%                'clk'       : 衛星時計推定結果グラフ
%                'clkb'      : 衛星時計推定結果表示(バイアス削除)
%                'clkf'      : 衛星時計推定結果一覧
%                'zpd'       : 対流圏天頂遅延推定結果グラフ
%                'zpdf'      : 対流圏天頂遅延推定結果一覧
%                'trg'       : 対流圏水平傾度推定結果グラフ
%                'rclk'      : 受信機時計グラフ
%                'rclkf'     : 受信機時計一覧
%                'rpos'      : 観測局位置結果一覧
%                'rpose'     : 観測局位置誤差一覧
%                'rposx'     : 観測局位置誤差XYZグラフ
%                'rpost'     : 観測局位置誤差ENUグラフ
%                'rposg'     : 観測局水平位置変位グラフ
%                'rposh'     : 観測局水平位置変位地図
%                'rposv'     : 観測局垂直位置変位地図
%                'bcp'       : 搬送波位相バイアス推定結果グラフ
%                'erp'       : 地球回転パラメータ推定結果グラフ
%                'stat'      : 統計情報
%                'res'       : 残差ENUグラフ
%                'resa'      : 残差AZELグラフ
%                'sats',sats : 衛星指定 {sat1,sat2,...}
%                'rcvs',rcvs : 観測局指定 {rcv1,rcv2,...}
%                'refclk',refclk : 時計基準衛星/観測局
%                'tspan',tspan : 時間範囲指定 [tstart,tend] (H)
%                'range',range : 値範囲指定 [vmin,vmax]
%                'dirs',dirs : ディレクトリ指定
%                'nsig',nsig : 標準偏差表示 n×σ
%                'jst'       : 時刻JST表示
%                'back'      : backward
%                'map',map   : 地図範囲 [lonw,lone,lats,latn]
% [argout] : なし
% [note]   :
% [version]: $Revision: 12 $ $Date: 2008-11-25 10:02:15 +0900 (轣ｫ, 25 11 2008) $
% [history]: 04/06/03  0.1  新規作成
%-------------------------------------------------------------------------------
prm=loadprm('prm_gpsest');
[td,ts]=caltomjd(prm.tstart); time=ts:prm.tint:ts+prm.tspan*3600;
sats=prm.sats; rcvs=prm.rcvs; refclk=''; types={}; 
dirs=prm.dirs; tspan=[0,inf]; range=[]; nsig=[]; utc_tai=prm.utc_tai; jst=0; n=1;
fb='f'; map=[137.3,140.7,36.4,38.6];
while n<=length(varargin)
    switch varargin{n}
    case 'td',    td    =varargin{n+1}; n=n+2;
    case 'time',  time  =varargin{n+1}; n=n+2;
    case 'sats',  sats  =varargin{n+1}; if ischar(sats), sats={sats}; end, n=n+2;
    case 'rcvs',  rcvs  =varargin{n+1}; if ischar(rcvs), rcvs={rcvs}; end, n=n+2;
    case 'refclk',refclk=varargin{n+1}; n=n+2;
    case 'tspan', tspan =varargin{n+1}; n=n+2;
    case 'range', range =varargin{n+1}; n=n+2;
    case 'dirs',  dirs  =varargin{n+1}; n=n+2;
    case 'nsig',  nsig  =varargin{n+1}; n=n+2;
    case 'jst',   jst=1; n=n+1;
    case 'back',  fb='b'; n=n+1;
    otherwise, types={types{:},varargin{n}}; n=n+1; end
end
if isempty(types), types={'eph'}; end
for n=1:length(types)
    if any(strcmp(types{n},{'eph','pos','vel','clk','clkb','srp','srpf'}))
        if isempty(refclk), refclk=prm.clkref; end
        for m=1:length(sats)
            PlotEstData(types{n},td,time,sats{m},{},refclk,tspan,range,utc_tai,nsig,jst,dirs,fb);
        end
    elseif any(strcmp(types{n},{'rclk','zpd','trg','rposx','rpost','rposg','bcp'}))
        if isempty(refclk), refclk=prm.clkref; end
        for m=1:length(rcvs)
            PlotEstData(types{n},td,time,{},rcvs{m},refclk,tspan,range,utc_tai,nsig,jst,dirs,fb);
        end
    elseif any(strcmp(types{n},{'rposh','rposv'}))
        mapdraw(map(1:2),map(3:4));
        lon=map(1)+(map(2)-map(1))*0.1;
        lat=map(3)+(map(4)-map(3))*0.9;
        if strcmp(types{n},'rposh')
            type='水平';
            mapplot(lon,lat,'arrow',0.05*3,0,'linewidth',2);
            maptext(lon+0.1,lat+0.05,'5cm');
        else
            type='垂直';
            mapplot(lon,lat-0.1,'arrow',0,0.05*3,'color','r','linewidth',2);
            maptext(lon+0.05,lat,'5cm');
        end
        if isinf(tspan(1)), tspan(1)=time(1)/3600; end
        if isinf(tspan(2)), tspan(2)=time(end)/3600; end
        for m=1:length(rcvs)
            PlotEstData(types{n},td,time,{},rcvs{m},refclk,tspan,range,utc_tai,nsig,jst,dirs,fb);
        end
        ts=tspan; if jst, ts=ts+9; st='JST'; else st='GPST'; end
        dt1=mjdtocal(td,ts(1)*3600); dt2=mjdtocal(td,ts(2)*3600);
        title(sprintf('観測点%s位置変位 : %04d/%02d/%02d %02d:%02d->%04d/%02d/%02d %02d:%02d %s',...
              type,dt1(1:5),dt2(1:5),st))
    elseif strcmp(types{n},'ephf')
        dt1=mjdtocal(td,time(1)+tspan(1)*3600); dt2=mjdtocal(td,time(end)-1);
        disp(sprintf('GPS ORBIT ESTIMATE ERROR : %04d/%02d/%02d %02d:%02d - %04d/%02d/%02d %02d:%02d',...
             dt1(1:5),dt2(1:5)))
        disp('       :         POSITION RMS ERR(m)                    VELOCITY RMS ERR(m/sec)')
        disp(' SAT   :   3DRMS   RADIAL ALONG-TRK CROSS-TRK   3DRMS     RADIAL    ALONK-TRK  CROSS-TRK')
        disp('----------------------------------------------------------------------------------------')
        rms=[];
        for m=1:length(sats)
            rm=PlotEstData(types{n},td,time,sats{m},{},[],tspan,range,utc_tai,nsig,jst,dirs,fb);
            if ~isempty(rm), rms=[rms;norm(rm(1:3)),rm(1:3),norm(rm(4:6)),rm(4:6)]; end
        end
        disp('----------------------------------------------------------------------------------------')
        i=find(~isnan(rms(:,1))); m1=mean(rms(i,:),1); m2=median(rms(i,:),1);
        disp(sprintf('AVERAGE: %7.4f %8.4f %8.4f %8.4f  %10.2E %10.2E %10.2E %10.2E',m1));
        disp(sprintf('MEDIAN : %7.4f %8.4f %8.4f %8.4f  %10.2E %10.2E %10.2E %10.2E',m2));
    elseif any(strcmp(types{n},{'clkf','rclkf','zpdf'}))
        if isempty(refclk), refclk=prm.clkref; end
        dt1=mjdtocal(td,time(1)+tspan(1)*3600); dt2=mjdtocal(td,time(end)-1);
        switch types{n}
        case 'clkf',  msg='SATELLITE CLOCK BIAS'; name='SAT';
        case 'rclkf', msg='RECEIVER CLOCK BIAS'; name='STA';
        case 'zpdf',  msg='TROPOSPHERIC ZENITH TOTAL DELAY'; name='STA';
        end
        disp(sprintf('%s ESTIMATION ERROR : %04d/%02d/%02d %02d:%02d - %04d/%02d/%02d %02d:%02d',...
             msg,dt1(1:5),dt2(1:5)))
        if strcmp(types{n},'clkf'), disp(sprintf('refclk = %s',refclk)), end
        disp('               ERROR(m)')
        disp(sprintf(' %s   :    BIAS      RMS',name))
        disp('-----------------------------')
        rms=[];
        if strcmp(types{n},'clkf')
            for m=1:length(sats)
                rm=PlotEstData(types{n},td,time,sats{m},{},refclk,tspan,range,utc_tai,nsig,jst,dirs,fb);
                if ~isempty(rm), rms=[rms;rm]; end
            end
        else
            for m=1:length(rcvs)
                rm=PlotEstData(types{n},td,time,{},rcvs{m},refclk,tspan,range,utc_tai,nsig,jst,dirs,fb);
                if ~isempty(rm), rms=[rms;rm]; end
            end
        end
        disp('-----------------------------')
        i=find(~isnan(rms(:,1)));
        disp(sprintf('AVERAGE: %8.4f %8.4f',mean(rms(i,:),1)))
        disp(sprintf('MEDIAN : %8.4f %8.4f',median(rms(i,:),1)))
    elseif any(strcmp(types{n},'rpose'))
        dt1=mjdtocal(td,time(1)+tspan(1)*3600); dt2=mjdtocal(td,time(end)-1);
        disp(sprintf('STATION POSITION ESTIMATION ERROR : %04d/%02d/%02d %02d:%02d - %04d/%02d/%02d %02d:%02d',...
             dt1(1:5),dt2(1:5)))
        disp(' STA   :   3D(m)      X(m)     Y(m)     Z(m)       E(m)     N(m)     U(m)')
        disp('-------------------------------------------------------------------------------------')
        rms=[];
        for m=1:length(rcvs)
            rm=PlotEstData(types{n},td,time,{},rcvs{m},[],tspan,range,utc_tai,nsig,jst,dirs);
            if ~isempty(rm), rms=[rms;rm]; end
        end
        disp('-------------------------------------------------------------------------------------')
        rms=rms(find(~isnan(rms(:,1))),:);
        disp(sprintf('RMS ERR:%8.4f   %8.4f %8.4f %8.4f   %8.4f %8.4f %8.4f',sqrt(mean(rms.^2,1))))
        disp(sprintf('AVERAGE:%8.4f   %8.4f %8.4f %8.4f   %8.4f %8.4f %8.4f',mean(rms,1)))
    elseif any(strcmp(types{n},'rpos'))
        dt1=mjdtocal(td,time(1)+tspan(1)*3600); dt2=mjdtocal(td,time(end)-1);
        disp(sprintf('STATION POSITION ESTIMATED : %04d/%02d/%02d %02d:%02d - %04d/%02d/%02d %02d:%02d',...
             dt1(1:5),dt2(1:5)))
        disp(' STA   :      X(m)           Y(m)           Z(m)       SDX(m)  SDY(m)  SDZ(m)')
        disp('-------------------------------------------------------------------------------------')
        rms=[];
        for m=1:length(rcvs)
            rm=PlotEstData(types{n},td,time,{},rcvs{m},[],tspan,range,utc_tai,nsig,dirs,fb);
        end
        disp('-------------------------------------------------------------------------------------')
    elseif any(strcmp(types{n},'erp'))
        PlotEstData(types{n},td,time,{},{},[],tspan,range,utc_tai,nsig,jst,dirs,fb)
    elseif any(strcmp(types{n},'stat'))
        PlotEstStat(td,time,sats,rcvs,tspan,dirs,fb)
    elseif any(strcmp(types{n},'res'))
        PlotResidual(td,time,sats,rcvs,tspan,dirs,fb)
    elseif any(strcmp(types{n},'resa'))
        PlotResidualAzel(td,time,sats,rcvs,tspan,dirs,fb)
    end
end

% 推定データ表示 ---------------------------------------------------------------
function rms=PlotEstData(type,td,time,sat,rcv,refclk,tspan,range,utc_tai,nsig,jst,dirs,fb)
rms=[];
if ~isempty(rcv)
    sta=prm_gpsrcvs;
    i=find(strcmp(rcv,sta(:,1)));
    if ~isempty(i), rname=sta{i,3}; else rname=''; end
end

if any(strcmp(type,{'eph','ephf','pos','vel'}))&~isempty(sat)
    [tde,timee,xs,covs,types,prm]=LoadEst(td,time,['eph',fb],sat,tspan,dirs);
    [tdr,timer,eph,clk,zpd]=LoadRef(td,time,sat,dirs);
    if isempty(timee), return, end
    timer=timer+(tde-tdr)*86400;
    [time,i,j]=intersect(timee,timer);
    if prm.est.erp==1
        [td,tt,erp]=LoadEst(td,time,'erp','',tspan,dirs);
        erp(:,1:2)=erp(:,1:2)*pi/180/3600;
        for k=i'
            tut=td+(timee(k)+19+utc_tai)/86400;
            erpv=[ErpVar(tut,utc_tai),0,0];
            erpr=ReadErp(tut,prm.dirs.erp,prm.src.erp)+erpv;
            erpe=[erp(k,:),erpr(4:5)]+erpv;
            U=ecsftoecef(tut,erpr,utc_tai)'*ecsftoecef(tut,erpe,utc_tai);
            xs(k,1:3)=(U*xs(k,1:3)')';
            xs(k,4:6)=(U*xs(k,4:6)')';
        end
    end
    xs=xs(i,1:6);
    err=xs-eph(j,:);
    if isinf(tspan(2)), tspan(2)=time(end)/3600; end
    
    if strcmp(type,'eph')
        figure, subplot(2,1,1), hold on, grid on
        errp=sqrt(sum(err(:,1:3).^2,2));
        title([sat,' Position/Velocity Estimation Error'])
        plot(time/3600,errp)
        if ~isempty(covs)
            sigp=sqrt(sum(covs(i,1:3),2));
            plot(time/3600,2*sigp,'r:')
        end
        if isempty(range), rangep=[0,0.3]; else rangep=range; end
        ylabel('pos error(m)'), TimeAxis(td,time,tspan,rangep,jst)
        DrawText(sprintf('RMS : %.5fm',RmsErr(errp)))
        
        subplot(2,1,2), hold on, grid on
        errv=sqrt(sum(err(:,4:6).^2,2));
        plot(time/3600,errv)
        if ~isempty(covs)
            sigv=sqrt(sum(covs(i,4:6),2));
            plot(time/3600,2*sigv,'r:')
        end
        if isempty(range), rangev=[0,1E-4]; else rangev=range; end
        if isempty(nsig), nsig=2; end
        ylabel('vel error(m)'), TimeAxis(td,time,tspan,rangev,jst)
        DrawText(sprintf('RMS : %.5Em/sec',RmsErr(errv)))
    
    elseif strcmp(type,'ephf')
        errp=zeros(length(time),3); errv=errp;
        for n=1:length(time)
            errp(n,:)=err(n,1:3)*ecsftosatf(xs(n,1:6)')';
            errv(n,:)=err(n,4:6)*ecsftosatf(xs(n,1:6)')';
        end
        errpv=[errp,errv]; errpv=errpv(find(~isnan(errpv(:,1))),:);
        if ~isempty(errpv)
            rms=sqrt(mean(errpv.^2,1));
            disp(sprintf('%-7s: %7.4f %8.4f %8.4f %8.4f  %10.2E %10.2E %10.2E %10.2E',...
                 sat,norm(rms(1:3)),rms(1:3),norm(rms(4:6)),rms(4:6)))
        else
            disp(sprintf('%-7s:     -       -         -       -          -          -          -          -     ',...
                 sat))
        end

    elseif strcmp(type,'pos')
        errp=zeros(length(time),3);
        for n=1:length(time), errp(n,:)=err(n,1:3)*ecsftosatf(xs(n,1:6)')'; end
        if isempty(range), range=[-0.5,0.5]; end
        figure, subplot(4,1,1), hold on, grid on
        errt=sqrt(sum(errp(:,1:3).^2,2));
        plot(time/3600,errt)
        if ~isempty(covs)
            sigp=sqrt(sum(covs(i,1:3),2));
            plot(time/3600,2*sigp,'r:')
        end
        ylabel('3drms(m)'), TimeAxis(td,time,tspan,[0,range(2)*2],jst)
        DrawText(sprintf('RMS : %.5fm',RmsErr(errt)))
        title([sat,' Position Estimation Error'])
        subplot(4,1,2), hold on, grid on, plot(time/3600,errp(:,1))
        ylabel('radial(m)'), TimeAxis(td,time,tspan,range,jst)
        DrawText(sprintf('RMS : %.5fm',RmsErr(errp(:,1))))
        subplot(4,1,3), hold on, grid on, plot(time/3600,errp(:,2))
        ylabel('along-track(m)'), TimeAxis(td,time,tspan,range,jst)
        DrawText(sprintf('RMS : %.5fm',RmsErr(errp(:,2))))
        subplot(4,1,4), hold on, grid on, plot(time/3600,errp(:,3))
        ylabel('cross-track(m)'), TimeAxis(td,time,tspan,range,jst)
        DrawText(sprintf('RMS : %.5fm',RmsErr(errp(:,3))))

    elseif strcmp(type,'vel')
        errv=zeros(length(time),3);
        for n=1:length(time), errv(n,:)=err(n,4:6)*ecsftosatf(xs(n,1:6)')'; end
        if isempty(range), range=[-1E-4,1E-4]; end
        figure, subplot(4,1,1), hold on, grid on
        errt=sqrt(sum(errv(:,1:3).^2,2));
        plot(time/3600,errt)
        if ~isempty(covs)
            sigp=sqrt(sum(covs(i,4:6),2));
            plot(time/3600,2*sigp,'r:')
        end
        ylabel('3drms(m/sec)'), TimeAxis(td,time,tspan,[0,range(2)],jst)
        DrawText(sprintf('RMS : %.5Em',RmsErr(errt)))
        title([sat,' 速度推定誤差'])
        subplot(4,1,2), hold on, grid on, plot(time/3600,errv(:,1))
        ylabel('radial(m/sec)'), TimeAxis(td,time,tspan,range,jst)
        DrawText(sprintf('RMS : %.5Em/sec',RmsErr(errv(:,1))))
        subplot(4,1,3), hold on, grid on, plot(time/3600,errv(:,2))
        ylabel('along-track(m/sec)'), TimeAxis(td,time,tspan,range,jst)
        DrawText(sprintf('RMS : %.5Em/sec',RmsErr(errv(:,2))))
        subplot(4,1,4), hold on, grid on, plot(time/3600,errv(:,3))
        ylabel('cross-track(m/sec)'), TimeAxis(td,time,tspan,range,jst)
        DrawText(sprintf('RMS : %.5Em/sec',RmsErr(errv(:,3))))
     end

elseif strcmp(type,'srp')&~isempty(sat)
    [td,time,xs,covs]=LoadEst(td,time,['eph',fb],sat,tspan,dirs);
    if isempty(time), return, end
    if isinf(tspan(2)), tspan(2)=time(end)/3600; end
    if isempty(range), range=[-inf,inf]; end
    if isempty(nsig), nsig=2; end
    figure, nx=size(xs,2)-6;
    for n=1:nx
        subplot(nx,1,n), hold on, grid on
        plot(time/3600,xs(:,6+n))
        if ~isempty(covs)&nsig>0
            sig=sqrt(covs(:,6+n));
            plot(time/3600,xs(:,6+n)+nsig*sig,'r:')
            plot(time/3600,xs(:,6+n)-nsig*sig,'r:')
        end
        ylabel(sprintf('srp prm(%d)',n)), TimeAxis(td,time,tspan,range,jst)
        if n==1, title([sat,' 太陽輻射圧パラメータ']), end
    end

elseif strcmp(type,'srpf')&~isempty(sat)
    [td,time,xs,covs,t,prm]=LoadEst(td,time(end-1),['eph',fb],sat,tspan,dirs);
    if isempty(time), return, end
    i=find(strcmp(prm.sats,sat));
    disp(sprintf('''%s'',[%2d,%7.1f,%5.1f,%12.6f,%10.6f,%10.6f,%10.6f]',...
         sat,prm.sat.srp(i,1:3),xs(end,7:10)))    

elseif any(strcmp(type,{'clk','clkb','clkf','rclk','rclkf'}))
    if ~isempty(sat), name=sat; else name=rcv; end
    C=299792458; f1=1.57542E9; f2=1.22760E9; lam1=C/f1; lam2=C/f2;
    cif=[f1^2;-f2^2]/(f1^2-f2^2);
    [tde,timee,xs,covs]=LoadEst(td,time,['clk',fb],name,tspan,dirs);
    if isempty(timee), return, end
    [tdr,timer,eph,clks,zpd]=LoadRef(td,time,name,dirs);
    if ~isempty(refclk), [tdr,timer,eph,clkr,zpd]=LoadRef(td,time,refclk,dirs);
    if isempty(nsig), nsig=2; end
    else clkr=[]; end
    timer=timer+(tde-tdr)*86400;
    if isinf(tspan(2)), tspan(2)=timee(end)/3600; end
    if isempty(range), range=[-0.5,0.5]; end
    [time,i,j]=intersect(timee,timer);
    if ~isempty(clks), clks=clks(j,1); end
    if ~isempty(clkr), clks=clks-clkr(j,1); end
    errc=xs(i,1)-clks;
    [rms,ave,st]=RmsErr(errc); rms=[ave,rms];
    if strcmp(type,'clkf')|strcmp(type,'rclkf')
        if isempty(sat), name=rcv; else name=sat; end
        disp(sprintf('%-7s: %8.4f %8.4f',name,rms))
        return
    end
    figure
    subplot(2,1,1), hold on, grid on
    title([name,' 時計推定誤差 (REFCLK : ',refclk,')'])
    plot(timee/3600,xs(:,1))
    ylabel('clock bias(m)'), TimeAxis(td,time,tspan,[-inf,inf],jst)
    h=plot(time/3600,clks,'m'); legend(h,'IGS')
    subplot(2,1,2), hold on, grid on
    plot(time/3600,errc,'-','markersize',1)
    if nsig>0
        sig=sqrt(covs(i,1));
        plot(time/3600, nsig*sig,'r:')
        plot(time/3600,-nsig*sig,'r:')
    end
    ylabel('clock bias error(m)'), TimeAxis(td,time,tspan,range,jst)
    DrawText(sprintf('MEAN : %.5fm RMS : %.5fm',rms))

elseif any(strcmp(type,{'zpd','zpdf'}))&~isempty(rcv)
    [tde,timee,xs,covs]=LoadEst(td,time,['zpd',fb],rcv,tspan,dirs);
    if isempty(timee), return, end
    timer=(floor((time(1)-3600)/7200):floor((time(end)+3600)/7200))*7200+3600;
    zpd=readtrop(td,timer,rcv,dirs.trop);
    if isinf(tspan(2)), tspan(2)=timee(end)/3600; end
    if isempty(range), range=[-0.1,0.1]; end
    i=find(tspan(1)*3600<=timee&timee<tspan(2)*3600);
    ave=mean(xs(i,1));
    if isnan(ave), range=[-inf,inf]; else range=range+ave; end
    if isempty(nsig), nsig=0; end
    err=repmat(nan,length(timer),1);
    for n=1:length(timer)
        i=find(timer(n)-3600<=timee&timee<timer(n)+3600); % 前後1H間平均で比較
        if ~isempty(i), err(n)=mean(xs(i,1))-zpd(n); end
    end
    [rms,ave,st]=RmsErr(err); rms=[ave,rms];
    if strcmp(type,'zpdf')
        disp(sprintf('%-7s: %8.4f %8.4f',rcv,rms))
        return
    end
    figure, hold on, grid on
    title([rcv,' ',rname,' 対流圏遅延推定値'])
    plot(timee/3600,xs(:,1))
    if ~isempty(covs)&nsig>0
        sig=sqrt(covs(:,1));
        plot(timee/3600,xs(:,1)+nsig*sig,'r:')
        plot(timee/3600,xs(:,1)-nsig*sig,'r:')
    end
    ylabel('zenith path delay(m)'), TimeAxis(td,time,tspan,range,jst)
    if isempty(zpd), return, end
    plot(timer/3600,zpd,'m');
    h=plot(timer/3600,zpd,'m.'); legend(h,'IGS')
    DrawText(sprintf('RMS : %.5fm',rms))

elseif any(strcmp(type,'trg'))&~isempty(rcv)
    [tde,timee,xs,covs]=LoadEst(td,time,['zpd',fb],rcv,tspan,dirs);
    [tdr,timer,eph,clks,zpd]=LoadRef(td,time,rcv,dirs);
    if isempty(timee), return, end
    timer=timer+(tde-tdr)*86400;
    if isinf(tspan(2)), tspan(2)=timee(end)/3600; end
    if isempty(range), range=[-0.03,0.03]; end
    if isempty(nsig), nsig=2; end
    labels={'trop-grad(NS)','trop-grad(EW)'};
    
    figure
    for n=1:2
        subplot(2,1,n), hold on, grid on
        if n==1, title([rcv,' ',rname,' 対流圏水平傾度推定値']), end
        plot(timee/3600,xs(:,n+1))
        if ~isempty(covs)&nsig>0
            sig=sqrt(covs(:,n+1));
            plot(timee/3600,xs(:,n+1)+nsig*sig,'r:')
            plot(timee/3600,xs(:,n+1)-nsig*sig,'r:')
        end
        ylabel(labels{n}), TimeAxis(td,time,tspan,range,jst)
    end

elseif any(strcmp(type,{'rpos','rpose','rposx','rpost','rposg','rposh','rposv'}))&~isempty(rcv)
    [tde,time,xs,covs]=LoadEst(td,time,['pos',fb],rcv,tspan,dirs);
    if isempty(time), return, end
    if isinf(tspan(2)), tspan(2)=time(end)/3600; end
    if isempty(range), range=[-0.1,0.1]; end
    if isempty(nsig), nsig=2; end
    if any(strcmp(type,{'rpose','rposx','rpost','rposg'}))
        pos=ReadPos(td,time(end),rcv,dirs.pos)';
        err=xs-repmat(pos',size(xs,1),1);
    end
    gpos=eceftogeod(xs(end,:)');
    [epos,E]=geodtoecef(gpos);
    
    switch type
    case 'rpos',
        disp(sprintf('%-6s :%14.4f %14.4f %14.4f  %7.4f %7.4f %7.4f : %s',rcv,xs(end,:),sqrt(covs(end,:)),rname))
    case 'rpose',
        errg=err(end,:)*E';
        rms=[norm(err(end,:)),err(end,:),errg];
        disp(sprintf('%-6s :%8.4f   %8.4f %8.4f %8.4f   %8.4f %8.4f %8.4f  : %s',rcv,rms,rname))
    case 'rposx',
        label={'pos x (m)','pos y (m)','pos z (m)'};
        figure
        for n=1:3
            subplot(3,1,n), hold on, grid on
            if n==1, title([rcv,' ',rname,' 観測点位置']), end
            plot(time/3600,err(:,n))
            if ~isempty(covs)&nsig>0
                sig=sqrt(covs(:,n));
                plot(time/3600,err(:,n)+nsig*sig,'r:')
                plot(time/3600,err(:,n)-nsig*sig,'r:')
            end
            ylabel(label{n}), TimeAxis(td,time,tspan,range,jst)
        end
    case 'rpost'
        for n=1:size(err,1), err(n,:)=err(n,:)*E'; end
        label={'pos east (m)','pos north (m)','pos up (m)'};
        figure
        for n=1:3
            subplot(3,1,n), hold on, grid on
            if n==1, title([rcv,' ',rname,' 観測点位置']), end
            plot(time/3600,err(:,n))
            ylabel(label{n}), TimeAxis(td,time,tspan,range,jst)
        end
    case 'rposg'
        if jst, time=time+9*3600; end
        for n=1:size(err,1), err(n,:)=err(n,:)*E'; end
        figure, hold on
        plot(err(:,1),err(:,2))
        i=find(mod(time,3600)==0);
        h1=plot(err(i,1),err(i,2),'.'); h2=[];
        for i=find(mod(time,86400)==0)'
            h2=plot(err(i,1),err(i,2),'r.');
            dt=mjdtocal(td,time(i));
            text(err(i,1),err(i,2)-0.002,sprintf('%02d/%02d',dt(2:3)),...
                 'horizontalalignment','center','verticalalignment','top')
        end
        grid on, axis equal, ylim(range), xlim(range*1.2), legend(h2,'00:00 pos')
        xlabel('W-E(m)'), ylabel('S-N(m)')
        ts=tspan; if jst, ts=ts+9; st='JST'; else st='GPST'; end
        dt1=mjdtocal(td,ts(1)*3600); dt2=mjdtocal(td,ts(2)*3600);
        title(sprintf('%s %s 水平位置 : %04d/%02d/%02d %02d:%02d-%04d/%02d/%02d %02d:%02d %s',...
              rcv,rname,dt1(1:5),dt2(1:5),st))
    case {'rposh','rposv'}
        i=find((tspan(1)-3)*3600<=time&time<(tspan(1)+3)*3600); % 前後3H平均
        j=find((tspan(2)-3)*3600<=time&time<(tspan(2)+3)*3600);
        if ~isempty(i)&~isempty(j)
            dpos=E*(mean(xs(j,:))-mean(xs(i,:)))';
        else dpos=[0;0;0]; end
        mapplot(gpos(2),gpos(1),'k.');
        maptext(gpos(2),gpos(1)-1E-2,rname,'horizontalalignment','center',...
                'verticalalignment','top');
        if strcmp(type,'rposh')
            mapplot(gpos(2),gpos(1),'arrow',dpos(1)*3,dpos(2)*3,'linewidth',2);
        else
            mapplot(gpos(2),gpos(1),'arrow',0,dpos(3)*3,'color','r','linewidth',2);
        end
    end

elseif any(strcmp(type,'bcp'))&~isempty(rcv)
    [td,time,xs,covs,types]=LoadEst(td,time,['bcp',fb],rcv,tspan,dirs);
    if isempty(time), return, end
    figure, hold on, grid on
    title([rcv,' ',rname,' 搬送波位相バイアス推定値'])
    if isinf(tspan(2)), tspan(2)=time(end)/3600; end
    if isempty(range), range=[-inf,inf]; end
    if isempty(nsig), nsig=2; end
    plot(time/3600,xs)
    if nsig>0
        sig=sqrt(covs);
        plot(time/3600, nsig*sig,'r:')
        plot(time/3600,-nsig*sig,'r:')
    end
    ylabel('phase bias(m)'), TimeAxis(td,time,tspan,range,jst)

elseif any(strcmp(type,'erp'))
    [td,timee,xs,covs,types]=LoadEst(td,time,['erp',fb],'',tspan,dirs);
    if isempty(timee), return, end
    timer=timee(1):900:timee(end);
    for n=1:length(timer)
        erp(n,:)=ReadErp(td+(timer(n)+19+utc_tai)/86400,dirs.erp);
    end
    figure, hold on, grid on
    if isinf(tspan(2)), tspan(2)=timee(end)/3600; end
    if isempty(range), range=[-inf,inf]; end
    if isempty(nsig), nsig=2; end
    sig=sqrt(covs); labels={'xp(")','yp(")','UT1-UTC(sec)'};
    for n=1:3
        if n<3, erp(:,n)=erp(:,n)*180*3600/pi; end
        subplot(3,1,n), hold on, grid on
        h=plot(timer/3600,erp(:,n),'m'); if n==1, legend(h,'IERS'), end
        plot(timee/3600,xs(:,n))
        if nsig>0
            plot(timee/3600,xs(:,n)+nsig*sig(:,n),'r:')
            plot(timee/3600,xs(:,n)-nsig*sig(:,n),'r:')
        end
        ylabel(labels{n}), TimeAxis(td,time,tspan,range,jst)
        if n==1, title([' 地球回転パラメータ推定値']), end
    end
end

% 統計情報表示 -----------------------------------------------------------------
function PlotEstStat(td,time,sats,rcvs,tspan,dirs)

[td,t,rsats,rrcvs,indexs,ress]=LoadStat(td,time,dirs);

disp('           obs count    prefit rms postfit rms ')
disp('-----------------------------------------------')
for n=1:length(sats)
    i=find(strcmp(sats{n},rsats));
    if ~isempty(i)
        j=find((indexs(:,1)==i|indexs(:,2)==i)&tspan(1)*3600<=t&t<tspan(2)*3600);
        if ~isempty(j)
            rmsp=sqrt(mean(ress(j,1).^2));
            rmsf=sqrt(mean(ress(j,2).^2));
            disp(sprintf('%-6s : %10d %10.4f %10.4f',sats{n},length(j),rmsp,rmsf));
        end
    end
end
disp(' ')
disp('           obs count    prefit rms postfit rms ')
disp('-----------------------------------------------')
for n=1:length(rcvs)
    i=find(strcmp(rcvs{n},rrcvs));
    if ~isempty(i)
        j=find((indexs(:,3)==i|indexs(:,4)==i)&tspan(1)*3600<=t&t<tspan(2)*3600);
        if ~isempty(j)
            rmsp=sqrt(mean(ress(j,1).^2));
            rmsf=sqrt(mean(ress(j,2).^2));
            disp(sprintf('%-6s : %10d %10.4f %10.4f',rcvs{n},length(j),rmsp,rmsf));
        end
    end
end

% 残差表示(ENU) ---------------------------------------------------------------
function PlotResidual(td,time,sats,rcvs,tspan,dirs)

ti=sprintf('%04d/%02d/%02d %02d:%02d:%02.0f-',mjdtocal(td));
[td,t,rsats,rrcvs,indexs,ress]=LoadStat(td,time,dirs);

for n=1:length(rcvs)
    i=find(strcmp(rcvs{n},rrcvs));
    if ~isempty(i)
        j=find(indexs(:,2)==i&tspan(1)*3600<=t&t<tspan(2)*3600&abs(ress)<1);
        if ~isempty(j)
            re=ress(j).*cos(azels(j,2)).*sin(azels(j,1));
            rn=ress(j).*cos(azels(j,2)).*cos(azels(j,1));
            ru=ress(j).*sin(azels(j,2));
            figure
            subplot(3,1,1), plot(t(j)/3600,re,'.','markersize',1), ylabel('east(m)')
            axis([tspan,-0.1,0.1]), grid on, title([rcvs{n},' 残差 : ',ti])
            DrawText(sprintf('BIAS : %.4fm RMS : %.4fm',mean(re),sqrt(mean(re.^2))))
            subplot(3,1,2), plot(t(j)/3600,rn,'.','markersize',1), ylabel('north(m)')
            axis([tspan,-0.1,0.1]), grid on,
            DrawText(sprintf('BIAS : %.4fm RMS : %.4fm',mean(rn),sqrt(mean(rn.^2))))
            subplot(3,1,3), plot(t(j)/3600,ru,'.','markersize',1), ylabel('up(m)')
            axis([tspan,-0.1,0.1]), grid on, xlabel('time(H)')
            DrawText(sprintf('BIAS : %.4fm RMS : %.4fm',mean(ru),sqrt(mean(ru.^2))))
        end
    end
end

% 残差表示(AZEL) ---------------------------------------------------------------
function PlotResidualAzel(td,time,sats,rcvs,tspan,dirs)

ti=sprintf('%04d/%02d/%02d %02d:%02d:%02.0f-',mjdtocal(td));
[td,t,rsats,rrcvs,indexs,ress]=LoadStat(td,time,dirs);

for n=1:length(rcvs)
    i=find(strcmp(rcvs{n},rrcvs));
    if ~isempty(i)
        j=find(indexs(:,2)==i&tspan(1)*3600<=t&t<tspan(2)*3600&abs(ress)<1);
        if ~isempty(j)
            figure
            subplot(2,1,1), plot(azels(j,1)*180/pi,ress(j),'.','markersize',1), ylabel('residual(m)')
            axis([-180,180,-0.1,0.1]), grid on, xlabel('azimuth(deg)')
            title([rcvs{n},' 残差 : ',ti])
            subplot(2,1,2), plot(azels(j,2)*180/pi,ress(j),'.','markersize',1), ylabel('residual(m)')
            axis([0,90,-0.1,0.1]), grid on, xlabel('elevation(deg)')
        end
    end
end


% 推定データ読み込み -----------------------------------------------------------
function [td,t,xs,vars,types,prm]=LoadEst(td,time,type,name,tspan,dirs)
t=[]; xs=[]; vars=[]; types={}; prm=[]; sats={}; tunit=86400;
ts=floor(time(1)/tunit)*tunit:tunit:time(end)-1;
for n=1:length(ts)
    dt=mjdtocal(td,ts(n));
    file=fullfile(dirs.est,sprintf('%s_%s_%04d%02d%02d%02d.mat',type,...
                  name,dt(1:4)));
    if exist(file)
        load(file), %disp(['load : ',file])
        [tdd,tss]=caltomjd(epoch);
        if isempty(td), td=tdd; t=time+tss;
        else t=[t;time+(tdd-td)*86400+tss]; end
        xs=[xs;data]; vars=[vars;covs];
    end
end
i=find(tspan(1)*3600<=t&t<tspan(2)*3600);
t=t(i); xs=xs(i,:);
if ~isempty(vars), vars=vars(i,:); end

% 参照データ読み込み -----------------------------------------------------------
function [tdd,t,ephr,clkr,zpdr]=LoadRef(tdd,time,name,dirs);
C=299792458; t=[]; ephr=[]; clkr=[]; zpdr=[]; tunit=10800;
ts=floor(time(1)/tunit)*tunit:tunit:time(end)-1;
for n=1:length(ts)
    dt=mjdtocal(tdd,ts(n));
    file=fullfile(dirs.ref,sprintf('ref_%04d%02d%02d%02d.mat',dt(1:4)));
    if exist(file)
        load(file)
        t=[t;time(:)+(td-tdd)*86400];
        n=find(strcmp(name,sats));
        m=find(strcmp(name,rcvs));
        if ~isempty(n)
            ephr=[ephr;eph(:,:,n)];
            clkr=[clkr;C*dts(:,n)];
        elseif ~isempty(m)
            clkr=[clkr;C*dtr(:,m)];
            zpdr=[zpdr;zpd(:,m)];
        end
    else, disp(['no file : ',file]), end
end

% 統計データ読み込み -----------------------------------------------------------
function [td,t,sats,rcvs,indexs,ress,azels]=LoadStat(td,time,dirs)
t=[]; sats={}; rcvs={}; indexs=[]; ress=[]; tunit=10800;
ts=floor(time(1)/tunit)*tunit:tunit:time(end)-1;
for n=1:length(ts)
    dt=mjdtocal(td,ts(n));
    file=fullfile(dirs.est,sprintf('res_%04d%02d%02d%02d.mat',dt(1:4)));
    if exist(file)
        load(file), %disp(['load : ',file])
        i=find(~outl);
        [tdd,tss]=CalToMjd(epoch);
        if isempty(td), td=tdd; t=time(i)+tss;
        else t=[t;time(i)+(tdd-td)*86400+tss]; end
        sats=prm.sats;
        rcvs=prm.rcvs;
        indexs=[indexs;index(i,:)];
        ress=[ress;residual(i,:)];
    end
end

% 時間軸表示 -------------------------------------------------------------------
function TimeAxis(td,time,tspan,range,jst)
if isinf(tspan(1)), tspan(1)=time(1)/3600; end
if isinf(tspan(2)), tspan(2)=time(end)/3600; end
ylim(range)
taxis(td,tspan,jst)

% テキスト描画 ---------------------------------------------------------------
function DrawText(string,pos)
if nargin<2, pos=2; end
p={'center','right'};
x=get(gca,'xlim'); y=get(get(gca,'title'),'position');
text(x(2),y(2),string,'horizontalalignment',p{pos})

% RMS誤差 --------------------------------------------------------------------
function [rms,ave,st]=RmsErr(x)
i=find(~isnan(x));
if isempty(i), rms=nan; ave=nan; st=nan; return, end
rms=sqrt(mean(x(i).^2));
ave=mean(x(i));
st=std(x(i));
if isempty(rms), rms=0; end

