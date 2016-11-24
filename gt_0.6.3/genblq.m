function genblq(file,td,ts,rcvs,dirs,src,path,tunit)
%-------------------------------------------------------------------------------
% [system] : GpsTools
% [module] : generate ocean loading parameters
% [func]   : generate ocean loading parameters
% [argin]  : file  = output file
%            td    = date (mjd)
%            ts    = time (sec)
%            rcvs  = stations list
%            dirs  = station position data directory
%            src   = station position source
%            path  = gotic2 path
%            tunit = station position processing unit time (hr)
% [argout] : none
% [note]   : BLQ format, 11 constituents
% [version]: $Revision: 8 $ $Date: 06/07/20 11:20 $
%            Copyright(c) 2004-2006 by T.Takasu, all rights reserved
% [history]: 04/09/19  0.1  new
%            04/12/11  0.2  support gui interface
%            06/07/20  0.3  add argin tunit
%-------------------------------------------------------------------------------
wave={'M2','S2','N2','K2','K1','O1','P1','Q1','Mf','Mm','Ssa'};
omodel='nao';

if nargin<5, dirs=''; end
if nargin<6, src =''; end
if nargin<7, path=''; end
if nargin<8, tunit=24; end
if ischar(rcvs), rcvs={rcvs}; end

fo=fopen(file,'wt'); if fo==-1, return, end

gut('newmsgbar','msgbar','',[360,50],[0,0.7,1],1);

WriteHead(fo,wave,omodel)

for n=1:length(rcvs)
    msg=['generating ocean loading parameters : ',rcvs{n}];
    if gut('setmsgbar','msgbar',msg,(n-1)/length(rcvs)), break; end
    rpos=readpos(td,ts,rcvs{n},dirs,src,tunit)';
    if ~isnan(rpos(1))
        gpos=eceftogeod(rpos);
        [odisp,ophas]=otidedisp(gpos,wave,omodel,path);
        WriteBody(fo,rcvs{n},gpos,wave,odisp,ophas)
    end
end
fclose(fo);
gut('closemsgbar','msgbar');

% write blq file header --------------------------------------------------------
function WriteHead(fo,wave,omodel)

fprintf(fo,'$$ Ocean loading displacement\n');
fprintf(fo,'$$\n');
fprintf(fo,'$$ Calculated by GOTIC2 (%s ver.%s)\n',mfilename,gpstools('version'));
fprintf(fo,'$$\n');
fprintf(fo,'$$ COLUMN ORDER:');
for n=1:length(wave), fprintf(fo,'%4s',wave{n}); end
fprintf(fo,'\n');
fprintf(fo,'$$\n');
fprintf(fo,'$$ Ocean tide model: %s\n',omodel);
fprintf(fo,'$$\n');
fprintf(fo,'$$ END HEADER\n');

% write blq file body ----------------------------------------------------------
function WriteBody(fo,rcv,gpos,wave,odisp,ophas)

fprintf(fo,'$$\n');
fprintf(fo,'  %s\n',rcv);
str=sprintf('%s,',rcv);
lon=gpos(2); if lon<0, lon=lon+360; end
fprintf(fo,'$$ %-36s RADI TANG lon/lat:  %8.4f  %8.4f\n',str,lon,gpos(1));
for n=1:3
    fprintf(fo,' ');
    for m=1:length(wave)
        str=sprintf('%7.5f',odisp(n,m));
        fprintf(fo,' %s',str(2:end));
    end
    fprintf(fo,'\n');
end
for n=1:3
    fprintf(fo,' ');
    for m=1:length(wave), fprintf(fo,'%7.1f',ophas(n,m)); end
    fprintf(fo,'\n');
end
