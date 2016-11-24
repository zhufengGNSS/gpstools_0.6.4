function clktornx(tstart,tend,tint,tunit,dirs,fb,comment)
%-------------------------------------------------------------------------------
% [system] : GpsTools
% [module] : convert clock estimation to rinex-clock format
% [func]   : convert clock estimation to rinex-clock format
% [argin]  : tstart : start time [year,month,day,hour,min,sec]
%            tend   : end time   [year,month,day,hour,min,sec]
%            tint   : time interval (sec)
%            tunit  : unit time (hr)
%            dirs   : directries
%                     dirs.est = estimation data
%                     dirs.clk = output directory
%            (fb)     : forward/backward estimation
%                     ('clkf':forward,'clkb':backward,'clkfb':smoothed)
%            (comment: comment)
%                comment{1}: reference satellite orbit
%                comment{2}: satellite antenna model
%                comment{3}: comment
% [argout] : none
% [note]   :
% [version]: $Revision: 12 $ $Date: 2008-11-25 10:02:15 +0900 (火, 25 11 2008) $
%            Copyright(c) 2004-2006 by T.Takasu, all rights reserved
% [history]: 06/05/09  0.1  new
%-------------------------------------------------------------------------------
if nargin<5, dirs.est=''; dirs.clk=''; end
if nargin<6, fb='clkf'; end
if nargin<7, comment={'','',''}; end
sats={}; for n=1:31, sats={sats{:},sprintf('GPS%02d',n)}; end
fname='clk';

[td,ts]=caltomjd(tstart);
[tn,te]=caltomjd(tend);
te=te+(tn-td)*86400; tu=tunit*3600;

for t=(floor(ts/tu):floor(te/tu))*tu
    time=t:tint:t+tu-tint;
    [clks,sigs,refc]=readclk(td,time,sats,dirs.est,fb,tunit);
    if ~isempty(clks)
        tdd=td+floor(t/86400); hr=floor(mod(t,86400)/3600);
        gpsd=td-44244; gpsw=floor(gpsd/7);
        file=sprintf('%s%04d%1d_%02d.clk',fname,gpsw,floor(gpsd-gpsw*7),hr);
        file=fullfile(dirs.clk,file);
        f=fopen(file,'wt');
        if f>0
            i=[];
            for n=1:length(sats)
                if any(~isnan(clks(:,1,n))), i=[i,n]; end
            end
            WriteHeader(f,td,time,sats(i),gpsw,gpsd,refc,comment);
            WriteBody(f,td,time,sats(i),clks(:,:,i),sigs(:,:,i),refc);
            fclose(f);
        end
    end
end

% write sp3 header -------------------------------------------------------------
function WriteHeader(f,td,time,sats,gpsw,gpsd,refc,comment)
pgm=['GPSTOOLS VER.',gpstools('version')];
comment1=['REFERENCE SATELLITE ORBIT      : ',comment{1}];
comment2=['SATELLITE ANTENNA PHASE CENTER : ',comment{2}];
comment3=comment{3};
agency=''; center=''; date=datestr(now);
fprintf(f,'     2.00           C                                       RINEX VERSION / TYPE\n');
fprintf(f,'%-20.20s%-20.20s%-20.20sPGM / RUN BY / DATE \n',pgm,agency,date(1:end-3));
fprintf(f,'%-60.60sCOMMENT             \n',comment1);
fprintf(f,'%-60.60sCOMMENT             \n',comment2);
fprintf(f,'%-60.60sCOMMENT             \n',comment3);
fprintf(f,'     2    AR    AS                                          # / TYPES OF DATA   \n');
fprintf(f,'%-20.20s                                        ANALYSIS CENTER     \n',center);
fprintf(f,'     1                                                      # OF CLK REF        \n');
fprintf(f,'%-4.4s                                                        ANALYSIS CLK REF    \n',refc);
fprintf(f,'    %2d                                                      # OF SOLN STA / TRF \n',1);
fprintf(f,'%-4.4s %-20.20s%11d %11d %11dSOLN STA NAME / NUM \n',refc,'',[0;0;0]);
fprintf(f,'    %2d                                                      # OF SOLN SATS      \n',length(sats));
for n=1:length(sats)
    fprintf(f,'G%-2.2s ',sats{n}(4:end));
    if mod(n,15)==0, fprintf(f,'PRN LIST            \n'); end
end
m=15-mod(length(sats),15);
if m<15, for n=1:m, fprintf(f,'    '); end, fprintf(f,'PRN LIST            \n'); end

fprintf(f,'                                                            END OF HEADER       \n');

% write clk body ---------------------------------------------------------------
function WriteBody(f,td,time,sats,clks,sigs,refc)
for n=1:length(time)
    if any(any(~isnan(clks(n,:,:))))
        fprintf(f,'AR %-4.4s %s %2d   %s\n',refc,tstr(td,time(n)),2,cstr(0,0));
    end
    for m=1:length(sats)
        sat=['G',sats{m}(4:end)];
        if ~isnan(clks(n,:,m))
            fprintf(f,'AS %-4.4s %s %2d   %s\n',sat,tstr(td,time(n)),...
                    2,cstr(clks(n,:,m),sigs(n,:,m)));
        end
    end
end

% clock string ------------------------------------------------------------------
function s=cstr(clk,sig), s=sprintf('%20.12E %20.12E',clk,sig); s([18,39])=[];

% time string ------------------------------------------------------------------
function s=tstr(td,t), s=sprintf('%04d %02d %02d %02d %02d %9.6f',mjdtocal(td,t));
