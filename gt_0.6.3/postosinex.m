function postosinex(td,span,tint,rcvs,indir,outdir,fb,tunit,trprm)
%-------------------------------------------------------------------------------
% [system] : GpsTools
% [module] : convert position estimation to sinex format
% [func]   : convert position estimation to sinex format
% [argin]  : tbd
% [argout] : none
% [note]   :
% [version]: $Revision: 4 $ $Date: 06/07/20 10:53 $
%            Copyright(c) 2004-2006 by T.Takasu, all rights reserved
% [history]: 05/06/13  0.1  new
%-------------------------------------------------------------------------------
if nargin<5, indir=''; end
if nargin<6, outdir=''; end
if nargin<7, fb='f'; end
if nargin<8, tunit=24; end
if nargin<9, trprm=''; end
ver=gpstools('version');
if fb=='f', time=(1:span)*86400-tint; else time=(0:span-1)*86400; end
for n=1:length(time)
    [poss,covs]=readpos(td,time(n),rcvs,indir,['pos',fb],tunit);
    [poss,tr]=TrCorr(td+time(n)/86400,poss,rcvs,trprm);
    SaveSinex(td,time(n),poss(1,:,:),covs(1,:,:),rcvs,indir,outdir,fb,tr,trprm,ver);
end

% save sinex file --------------------------------------------------------------
function SaveSinex(td,time,poss,covs,rcvs,indir,outdir,fb,tr,trprm,ver)
formatver=2.1; agency='';
dt=mjdtocal(td,time(1));
file=sprintf('pos%s%04d%02d%02d.snx',fb,dt(1:3));
file=gfilepath(outdir,file,dt);
f=fopen(file,'wt');
if f<0, disp(['warning : sinex file open error : ',file]), return, end
fprintf(f,'%%=SNX %4.2f %-3.3s %s %-3.3s %s %s P %05d 2 S           \n',...
        formatver,agency,sinext(caltomjd(datevec(now))),agency,...
        sinext(td+time(1)/86400),sinext(td+time(end)/86400),length(rcvs));
fprintf(f,'*-------------------------------------------------------------------------------\n');
fprintf(f,'+FILE/RERERENCE\n');
fprintf(f,' %-18.18s %-60.60s\n','DESCRIPTION','STATION POSITION ESTIMATIONS');
fprintf(f,' %-18.18s %-60.60s\n','SOFTWARE',['GPSTOOLS VER.',ver]);
fprintf(f,' %-18.18s %-60.60s\n','INPUT',indir);
fprintf(f,' %-18.18s %-60.60s\n','INPUT',['pos',fb]);
fprintf(f,' %-18.18s %-60.60s\n','INPUT',trprm);
fprintf(f,'-FILE/REFERENCE\n');
fprintf(f,'*-------------------------------------------------------------------------------\n');
fprintf(f,'+SITE/ID\n');
fprintf(f,'*Code Pt __Domes__ T _Station Description__ _Longitude_ _Latitude__ _Height\n');
rprm=prm_gpsrcvs;
rpos=readpos(td,time(1),rcvs,'','approx');
for n=1:length(rcvs)
    rp=rprm(strcmp(rprm(:,1),rcvs{n}),:);
    gpos=eceftogeod(rpos(1,:,n)');
    if ~isempty(rp)&~isempty(gpos)
        if gpos(2)<0, gpos(2)=gpos(2)+360; end
        fprintf(f,' %-4.4s  A %-9.9s P %-22.22s %3.0f %2.0f %4.1f %3.0f %2.0f %4.1f %7.1f\n',...
            rcvs{n}(end-3:end),rcvs{n},rp{3},degtodms(gpos(2)),degtodms(gpos(1)),gpos(3));
    end
end
fprintf(f,'-SITE/ID\n');
fprintf(f,'*-------------------------------------------------------------------------------\n');
fprintf(f,'+SOLUTION/ESTIMATE\n');
fprintf(f,'*Index _Type_ Code Pt Soln _Ref_Epoch__ Unit S __Estimated Value____ _Std_Dev___\n');
prm={'STAX','STAY','STAZ'}; index=1;
for n=1:length(rcvs)
    if all(~isnan(poss(1,:,n)))
        for m=1:3
            fprintf(f,' %5d %-6.6s %-4.4s  A ---- %s %-4.4s 2 %s %s\n',...
                    index,prm{m},rcvs{n}(end-3:end),sinext(td+time(1)/86400),'m',...
                    prmstr(poss(1,m,n),22),prmstr(sqrt(covs(1,m,n)),12));
            index=index+1;
        end
    end
end
fprintf(f,'-SOLUTION/ESTIMATE\n');
if ~isempty(tr)
    fprintf(f,'*-------------------------------------------------------------------------------\n');
    fprintf(f,'+SOLUTION/APRIORI\n');
    prm={'TX','TY','TZ','RX','RY','RZ','SC'}; units={'m','m','m','mas','mas','mas','ppb'};
    for n=1:7
        fprintf(f,' %5d %-6.6s %-4.4s  A ---- %s %-4.4s 2 %s %s\n',...
                n,prm{n},'----',sinext(floor(td+time(1)/86400)),units{n},...
                prmstr(tr(n),22),prmstr(0,12));
    end
    fprintf(f,'-SOLUTION/APRIORI\n');
end
fprintf(f,'*-------------------------------------------------------------------------------\n');
fprintf(f,'+SOLUTION/MATRIX_ESTIMATE\n');
fprintf(f,'-SOLUTION/MATRIX_ESTIMATE\n');
fprintf(f,'*-------------------------------------------------------------------------------\n');
fprintf(f,'%%ENDSNX\n');
fclose(f);
disp(['save : ',file])

% sinex time format ------------------------------------------------------------
function str=sinext(t)
dt=mjdtocal(t);
str=sprintf('%02d:%03d:%05d',mod(dt(1),100),floor(t)-caltomjd([dt(1),1,1])-1,...
            floor(mod(t,1)*86400));

% parameter string -------------------------------------------------------------
function s=prmstr(val,col)
if val==0, s=sprintf(['%.',num2str(col-5),'fE%+d'],0,0); s(1)=[];
else
    e=floor(log10(abs(val))+1);
    if val<0, s=sprintf(['%.',num2str(col-6),'fE%+d'],val/10^e,e); s(2)=[];
    else s=sprintf(['%.',num2str(col-5),'fE%+d'],val/10^e,e); s(1)=[]; end
end

% deg -> deg,min,sec -----------------------------------------------------------
function d=degtodms(deg), d=floor([deg,mod(deg*60,60),mod(deg*3600,60)]);

% tranform correction ----------------------------------------------------------
function [poss,tr]=TrCorr(t,poss,rcvs,trprm)
prm=[]; tr=[];
[d,f,e]=fileparts(trprm);
if strcmp(e,'.m'), wd=cd; cd(d); prm=feval(f); cd(wd); end
for n=1:size(prm,1)
    tt=caltomjd(prm(n,1:3)); if tt<=t&t<tt+1, tr=prm(n,4:end)'; break; end
end
if ~isempty(tr)
    for n=1:length(rcvs), poss(1,:,n)=helmert(poss(1,:,n),tr); end
end
