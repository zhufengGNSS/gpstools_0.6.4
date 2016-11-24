function [odisp,ophas]=otidedisp(gpos,wave,omodel,path)
%-------------------------------------------------------------------------------
% [system] : GpsTools
% [module] : generate ocean loading parameters
% [func]   : generate ocean loading parameters
% [argin]  : gpos = station position [lat,lon,h] (deg,m)
%            (wave)   = constituents {wv1,wv2,...}
%            (omodel) = tide model ('nao'=NAO.99b,'got'=GOTO,'csr'=CSR4.0)
%            (path)   = gotic2 path
% [argout] : odisp    = ocean loading amplitude (m)
%            ophas    = ocean loading phase (deg)
% [note]   : gotic2 front-end
%            gotic2 must be installded in argin path or <command path>\..\..\gotic2
%            GOTIC2 version 2001.05.16
%            reference:
%            Matsumoto, K., T. Sato, T. Takanezawa, and M. Ooe,
%            GOTIC2: A Program for Computation of Oceanic Tidal Loading Effect, 
%            J. Geod. Soc. Japan, 47, 243-248, 2001
% [version]: $Revision: 2 $ $Date: 06/07/08 14:12 $
%            Copyright(c) 2004-2006 by T.Takasu, all rights reserved
% [history]: 04/09/19  0.1  new
%            04/11/02  0.2  support matlab7
%-------------------------------------------------------------------------------
if nargin<2, wave={'M2','S2','N2','K2','K1','O1','P1','Q1','Mf','Mm','Ssa'}; end
if nargin<3, omodel='nao'; end
if nargin<4|isempty(path)
    [path,name]=fileparts(which(mfilename));
    path=[path,filesep,'..',filesep,'..',filesep,'gotic2']; % path of gotic2
end
fctr='tmp.ctr'; fout='tmp.out';
odisp=zeros(3,length(wave)); ophas=zeros(3,length(wave));

if exist(path,'dir')&exist(fullfile(path,'gotic2.exe'),'file')
    wd=pwd;
    cd(path);
    if GenCtr(fctr,gpos,wave,omodel)
        dos(['gotic2.exe <',fctr,' >',fout]);
        [odisp,ophas]=ReadOut(fout,wave);
        delete(fout);
    end
    delete(fctr);
    cd(wd);
else
    disp(['gotic2 not installed : ',path])
end

% output gotic2 control file ---------------------------------------------------
function stat=GenCtr(file,gpos,wave,omodel)
fo=fopen(file,'wt');
if fo==-1, disp(['control file open error : ',file]), stat=0; return, end
fprintf(fo,'STAPOSD 0000, %.5f, %.5f, %.3f, 0.0\n',gpos([2,1,3]));
wv=wave{1}; for n=2:length(wave), wv=[wv,',',wave{n}]; end
fprintf(fo,'WAVE %s\n',wv);       % constituents
fprintf(fo,'KIND RD,HD\n');       % output type
fprintf(fo,'GREENF 1\n');         % green function : 1=G-B,2=1066A
fprintf(fo,'POINTL ON\n');        % point mass calculation : ON,OFF
fprintf(fo,'MASSCON ON\n');       % mass conservation : ON,OFF
fprintf(fo,'OMODEL %s\n',omodel); % tide model
fprintf(fo,'GTIME 2004, 1, 1\n');
fprintf(fo,'PATAN2 OFF\n');
fprintf(fo,'VERBOUS OFF\n');
fprintf(fo,'END\n');
fclose(fo);
stat=1;

% read gotic2 output file ------------------------------------------------------
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
