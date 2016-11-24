% éûåvïœìÆâêÕ
function sample2

dirs.clk='d:\gps\clk';
dirs.pos='d:\gps\pos';

td=CalToMjd([2003,12,1]);
time=0:300:86400*31;
tau=[3E2,9E2,3E3,9E3,3E4,9E4];

sats={...
'GPS01','GPS02','GPS03','GPS04','GPS05','GPS06','GPS07','GPS08','GPS09','GPS10',...
'GPS11','GPS13','GPS14','GPS15','GPS16','GPS17','GPS18','GPS20','GPS21','GPS23',...
'GPS24','GPS25','GPS26','GPS27','GPS28','GPS29','GPS30','GPS31'};
rcvs={...
'ALBH','ALGO','ALIC','ALRT','AMC2','AOML','AREQ','ARTU','ASC1','AUCK','BAHR','BILI',...
'BOGT','BRAZ','BRMU','CAS1','CEDU','CHAT','CHUR','COCO','CORD','CRO1','DAKA','DARW',...
'DRAO','DUBO','DWH1','FAIR','FORT','GLPS','GODE','GOLD','GOUG','GRAZ','GUAM','GUAT',...
'HERS','HOB2','HRAO','IISC','INVK','IRKT','JPLM','KARR','KERG','KIT3','KOKB','KOSG',...
'KOUR','KSTU','LAE1','LHAS','LPGS','MAC1','MADR','MAS1','MATE','MCM4','MDO1','MEDI',...
'METS','MIZU','MKEA','NKLG','NLIB','NNOR','NOT1','NOUM','NRC1','NRIL','NTUS','NYAL',...
'OBET','OHI2','ONSA','PDEL','PERT','PETP','PIE1','PIMO','POTS','PRDS','QAQ1','REYK',...
'RIOG','SANT','SCH2','SFER','STJO','SUTH','SUWN','SYOG','THTI','THU3','TIDB','TIXI',...
'TLSE','TOW2','TRO1','TSKB','TWTF','UNSA','URUM','USUD','VILL','WES2','WSRT','WTZR',...
'WUHN','YAKT','YAR2','YELL','YSSK','ZAMB','ZECK'};

dt0=MjdToCal(td+time(1)/86400); dt1=MjdToCal(td+time(end)/86400);
disp('     GPS SATELLITE/IGS STATION CLOCK ALLAN-DEVIATION : SIGMA(tau) (sec/sec)')
disp(sprintf('%-8s: %04d/%02d/%02d %02d:%02d - %04d/%02d/%02d %02d:%02d','TIME',dt0(1:5),dt1(1:5)))
disp(sprintf('%-8s: %.2E %.2E %.2E %.2E %.2E %.2E','TAU(sec)',tau))
disp('------------------------------------------------------------------------------')

clks=ReadClk(td,time,sats,dirs.clk);
for n=1:length(sats)
    sig=repmat(nan,length(tau),1);
    for m=1:length(tau)
        clk=clks(1:tau(m)/300:end,:,n);
        dclk=(clk(1:end-1)-clk(2:end))/tau(m);
        i=find(~isnan(dclk));
        if~isempty(i), sig(m)=std(dclk(i)); end
    end
    disp(sprintf('%-8s: %.2E %.2E %.2E %.2E %.2E %.2E : %s',sats{n},sig,''))
end

clks=ReadClk(td,time,rcvs,'clk');
for n=1:length(rcvs)
    [apc,ecc,ant,rec]=ReadRcv(td,0,rcvs{n},dirs.pos);
    sig=repmat(nan,length(tau),1);
    for m=1:length(tau)
        clk=clks(1:tau(m)/300:end,:,n);
        dclk=(clk(1:end-1)-clk(2:end))/tau(m);
        i=find(~isnan(dclk));
        if~isempty(i), sig(m)=std(dclk(i)); end
    end
    disp(sprintf('%-8s: %.2E %.2E %.2E %.2E %.2E %.2E : %s',rcvs{n},sig,rec))
end
