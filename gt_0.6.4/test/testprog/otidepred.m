function [dpos,dposf]=OtidePred(td,time,gpos,tide,wave,omodel)
%-------------------------------------------------------------------------------
% [system] : GpsTools
% [module] : ocean loading prediction
% [func]   : predict site displacement by ocean loading
% [argin]  : td,time = date(mjd),time vector(sec) (utc)
%            gpos = station position [lat,lon,h] (deg,m)
%            (tide)   = 1:solid earth+ocean loading,2:solid only,3:ocean only
%            (wave)   = constituents {wv1,wv2,...}
%            (omodel) = tide model('nao'=NAO.99b,'got'=GOTO,'csr'=CSR4.0)
% [argout] : dpos  = site displacement(ecef)(m) [x,y,z;...]
%            dposf = site displacement(local tangental)(m) [east,north,up;...]
% [note]   : gotic2 front-end
%            gotic2 must be installed in <command path>\..\..\gotic2
%            GOTIC2 version 2001.05.16
%            reference:
%            Matsumoto, K., T. Sato, T. Takanezawa, and M. Ooe,
%            GOTIC2: A Program for Computation of Oceanic Tidal Loading Effect, 
%            J. Geod. Soc. Japan, 47, 243-248, 2001
% [version]: $Revision: 12 $ $Date: 2008-11-25 10:02:15 +0900 (火, 25 11 2008) $
% [history]: 04/10/08  0.1  new
%-------------------------------------------------------------------------------
if nargin<4, tide=1; end
if nargin<5, wave={'ALL'}; end
if nargin<6, omodel='nao'; end

fctr='tmp.ctr';
fout='tmp.out';
[path,name]=fileparts(which(mfilename));
dirs=[path,filesep,'..',filesep,'..',filesep,'gotic2']; % path of GOTIC2

if exist(dirs,'dir')&exist(fullfile(dirs,'gotic2.exe'),'file')
    wd=pwd;
    cd(dirs);
    if GenCtr(fctr,gpos,td,time,tide,'RD',fout,wave,omodel)
        [s,w]=dos(['gotic2.exe <',fctr]);
        dposr=dlmread(fout);
    end
    if GenCtr(fctr,gpos,td,time,tide,'HD',fout,wave,omodel)
        [s,w]=dos(['gotic2.exe <',fctr]);
        dposh=dlmread(fout);
    end
    delete(fout);
    delete(fctr);
    cd(wd);
else
    disp(['gotic2 not installed : ',dirs])
end
[epos,E]=geodtoecef(gpos);
dposf=[dposh(:,[3,2]),dposr(:,2)];
for n=1:length(time), dpos(n,:)=dposf(n,:)*E; end

% output GOTIC2 control file ---------------------------------------------------
function stat=GenCtr(file,gpos,td,time,tide,type,outfile,wave,omodel)
ts=mjdtocal(td,time(1));
te=mjdtocal(td,time(end));
ti=(time(2)-time(1))/60;
fo=fopen(file,'wt');
if fo==-1, disp(['control file open error : ',file]), stat=0; return, end
fprintf(fo,'STAPOSD 0000, %.5f, %.5f, %.3f, 0.0\n',gpos([2,1,3]));
wv=wave{1}; for n=2:length(wave), wv=[wv,',',wave{n}]; end
fprintf(fo,'WAVE %s\n',wv);       % constituents
fprintf(fo,'KIND %s\n',type);     % output type
%fprintf(fo,'GREENF 1\n');        % green function : 1=G-B,2=1066A
%fprintf(fo,'POINTL ON\n');       % point mass caliculation : ON,OFF
%fprintf(fo,'MASSCON ON\n');
%fprintf(fo,'OMODEL %s\n',omodel); % tide model
%fprintf(fo,'GTIME 2004, 1, 1\n');
%fprintf(fo,'PATAN2 OFF\n');
fprintf(fo,'PREDICT %d,%04d%02d%02d%02d%02d,%04d%02d%02d%02d%02d,%.0f\n',tide,...
        ts(1:5),te(1:5),ti);
fprintf(fo,'PREFMT 4,1\n');
fprintf(fo,'PREOUT %s\n',outfile);
fprintf(fo,'END\n');
fclose(fo);
stat=1;

% read GOTIC2 output file ------------------------------------------------------
function [odisp,ophas]=ReadOut(file,wave)
odisp=zeros(3,length(wave)); ophas=zeros(3,length(wave));
sec=0; n=0;
fi=fopen(file,'rt');
if fi==-1, disp(['output file open error : ',file]), return, end
while 1
    str=fgetl(fi); if ~isstr(str), break, end
    
    if findstr(str,'Theoretical Earth tide'), sec=0; n=0;
    elseif findstr(str,'Radial displacement'), sec=1;
    elseif findstr(str,'Tangential displacement'), sec=2;
    elseif sec==1&findstr(str,'Oceanic tidal effects'), n=1;
    elseif sec==2&findstr(str,'Oceanic tidal effects'), sec=3;
    elseif sec==3&findstr(str,'E-W') n=2;
    elseif sec==3&findstr(str,'N-S') n=3;
    elseif sec==3&findstr(str,'Azimuth') n=0;
    elseif 1<=n&n<=3
        m=find(strcmp(sscanf(str,'%s*%f'),wave));
        if ~isempty(m)
            v=sscanf(str(8:end),'%f');
            odisp(n,m)=v(9);
            ophas(n,m)=v(10);
            if 2<=n&n<=3, ophas(n,m)=ophas(n,m)+180; end % reverse phase(NS,EW)
            if 180<ophas(n,m), ophas(n,m)=ophas(n,m)-360; end
        end
    end
end
fclose(fi);
