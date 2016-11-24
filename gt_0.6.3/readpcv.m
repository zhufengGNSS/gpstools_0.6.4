function [apc1,apc2,apv1,apv2,stat]=readpcv(rcv,ant,file,time)
%-------------------------------------------------------------------------------
% [system] : GpsTools
% [module] : read satellite/receiver antenna parameters
% [func]   : read satellite/receiver antenna phase center offset and variations
% [argin]  : rcv    = station name
%            ant    = receiver antenna model name
%            (file) = antenna pcv parameter file path (default:'igs_01.pcv')
%            (time) = date-time (mjd)
% [argout] : apc1 = L1 phase center offset (m) [up,north,east]
%            apc2 = L2 phase center offset (m) [up,north,east]
%            apv1 = L1 phase center variation (m) (elev=0:5:90deg)
%            apv2 = L2 phase center variation (m) (elev=0:5:90deg)
%            stat = status (1:ok,0:no file/no antenna type)
% [note]   : antenna alias name file (prm_gpsants.m in search path) returning :
%                alias{n,1} = antenna alias name
%                alias{n,2} = antenna name in antenna pcv parameter file
% [version]: $Revision: 15 $ $Date: 06/07/13 8:05 $
%            Copyright(c) 2004-2006 by T.Takasu, all rights reserved
% [history]: 04/11/02   0.1  new
%            04/12/08   0.2  change argin pcvfile
%            05/05/31   0.3  support ngs pcvfile
%            05/07/02   0.4  return nan if no pcv found
%            05/09/19   0.5  support antenna alias name
%            06/03/23   0.6  support m-file antenna pcv parameters
%            06/04/21   0.7  support antex antenna exchange format
%                            support satellite antenna/add argin time
%-------------------------------------------------------------------------------
persistent file_ rcvs ants apcs1 apcs2 apvs1 apvs2 valid
if nargin<3, file='igs_01.pcv'; end
if nargin<4, time=[]; end

apc1=[0,0,0]; apc2=[0,0,0]; apv1=zeros(1,19); apv2=zeros(1,19);

% read antenna pcv paremter file
if ~strcmp(file,file_)
    rcvs={}; ants={};
    [ds,fs,ext]=fileparts(file);
    switch ext
    case '.m'
        wd=pwd; cd(ds); apcs=feval(fs); cd(wd);
        for n=1:size(apcs,1)
            rcvs{n}=apcs{n,1};
            apcs1(n,:)=[apcs{n,[3,4,2]}];
            apcs2(n,:)=[apcs{n,[6,7,5]}];
            apvs1(n,:)=zeros(1,19);
            apvs2(n,:)=zeros(1,19);
        end
    case '.atx'
        [ants,apcs1,apcs2,apvs1,apvs2,valid]=readantex(file);
    otherwise
        [ants,apcs1,apcs2,apvs1,apvs2]=ReadPcvFile(file);
        valid=[];
    end
    file_=file;
end
% antenna alias name
[dirs,f]=fileparts(which(mfilename));
alias=fullfile(dirs,'data','ants_alias.txt');
if exist(alias)
    [anta,antas]=textread(alias,'%s%s%*[^\n]','delimiter',',',...
                         'commentstyle','matlab');
    i=find(strcmp(deblank(ant),anta(:,1)));
    if ~isempty(i), ant=antas{i}; end
end
% search station name
i=find(strcmp(rcv,rcvs));

% search antenna type
if isempty(i)
    for n=1:length(ants)
        if isempty(valid)|isempty(time)|(valid(n,1)<=time&time<=valid(n,2))
            if findstr(ant,ants{n})==1, i=n; break; end
        end
    end
end
% search antenna without radome
if isempty(i)
    antr=strtok(ant);
    for n=1:length(ants)
        if isempty(valid)|isempty(time)|(valid(n,1)<=time&time<=valid(n,2))
            if findstr(antr,ants{n})==1, i=n; break; end
        end
    end
end
if isempty(i)
    stat=0;
elseif findstr(ant,'GPS')
    apc1=apcs1(i,:);
    apc2=apcs2(i,:);
    apv1=-apvs1(i,:);
    apv2=-apvs2(i,:);
    stat=1;
else
    apc1=apcs1(i,[3,1,2]);  % N/E/U->U/N/E
    apc2=apcs2(i,[3,1,2]);
    apv1=apvs1(i,end:-1:1); % 90:5:0->0:5:90
    apv2=apvs2(i,end:-1:1);
    stat=1;
end

% read antenna pcv parameter file  ---------------------------------------------
function [ants,apcs1,apcs2,apvs1,apvs2]=ReadPcvFile(file)
ants={}; apcs1=[]; apcs2=[]; apvs1=[]; apvs2=[];
f=fopen(file,'rt'); if f<0, return, end
n=0; m=0; p=1; s=16;
while 1
    str=fgetl(f); if ~isstr(str), break, end
    
    if m<=0|6<m
        p1=findstr(str,'MODEL #');
        p2=findstr(str,'ANTENNA ID');
        if     ~isempty(p1), p=p1; s=20;
        elseif ~isempty(p2), p=p2; s=16;
        elseif isletter(SubStr(str,1,1))
            n=n+1; m=1;
            ants{n,1}=deblank(SubStr(str,p,s));
            apcs1(n,1:3)=0; apcs2(n,1:3)=0; apcv1(n,1:19)=0; apcv2(n,1:19)=0;
        end
    else
        switch m
        case 1, for k= 1: 3, apcs1(n,k)=StrToNum(str,10*(k- 1)+1,10)*1E-3; end
        case 2, for k= 1:10, apvs1(n,k)=StrToNum(str, 6*(k- 1)+1, 6)*1E-3; end
        case 3, for k=11:19, apvs1(n,k)=StrToNum(str, 6*(k-11)+1, 6)*1E-3; end
        case 4, for k= 1: 3, apcs2(n,k)=StrToNum(str,10*(k- 1)+1,10)*1E-3; end
        case 5, for k= 1:10, apvs2(n,k)=StrToNum(str, 6*(k- 1)+1, 6)*1E-3; end
        case 6, for k=11:19, apvs2(n,k)=StrToNum(str, 6*(k-11)+1, 6)*1E-3; end
        end
        m=m+1;
    end
end
fclose(f);

% substring --------------------------------------------------------------------
function sstr=SubStr(str,pos,ns), sstr=str(pos:min(pos+ns-1,length(str)));

% string->number ---------------------------------------------------------------
function num=StrToNum(str,pos,ns)
num=sscanf(str(pos:min(pos+ns-1,length(str))),'%f');
if isempty(num), num=0; end
