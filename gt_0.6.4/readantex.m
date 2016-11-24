function [ants,apc1,apc2,apv1,apv2,valid]=readantex(file)
%-------------------------------------------------------------------------------
% [system] : GpsTools
% [module] : read antex antenna parameters file
% [func]   : read antex antenna parameters file
% [argin]  : file   = file path
% [argout] : ants   = receiver antenna type/satellite name [name,radome]
%            apc1   = L1 receiver anttena phase center offset  [n,e,u] (m)
%                     L1 satellite anttena phase center offset [x,y,z] (m)
%            apc2   = L2 receiver anttena phase center offset  [n,e,u] (m)
%                     L2 satellite anttena phase center offset [x,y,z] (m)
%            apv1   = L1 receiver/satellite anttena pcv (m) (73x19)
%            apv2   = L2 receiver/satellite anttena pcv (m) (73x19)
%                     apv?(i,j,k)=ants{i},az(j),ze(k)/nadir(k) pcv (m)
%                          az=0:5:360deg,ze=0:5:90deg/nadir=0:1:18deg
%            valid  = valid perod [start,end] (mjd)
% [note]   :
% [version]: $Revision: 16 $ $Date: 2008-12-12 15:49:30 +0900 (金, 12 12 2008) $
%            Copyright(c) 2004-2008 by T.Takasu, all rights reserved
% [history]: 06/04/22   0.1  new
%            08/11/24   0.2  change size of argout apv1,apv2
%                            support satellite GLONASS
%                            fix bug on error-stop if format error
%-------------------------------------------------------------------------------
ants={}; apc1=[]; apc2=[]; apv1=[]; apv2=[]; valid=[];

fd=fopen(file,'rt'); if fd<0, return, end

% read antex header
type=ReadHead(fd);

% read antex body
[ants,apc1,apc2,apv1,apv2,valid]=ReadData(fd);

fclose(fd);

% read antex header ------------------------------------------------------------
function type=ReadHead(fd)
type='';
while 1
    str=fgetl(fd); if ~isstr(str), break; end
    label=SubStr(str,61,20);
    if findstr(label,'ANTEX VERSION / SYST')
        type=SubStr(str,21,1);
    elseif findstr(label,'END OF HEADER'), break, end
end

% read rinex clock body --------------------------------------------------------
function [ants,apc1,apc2,apv1,apv2,valid]=ReadData(fd)
ants={}; apc1=[]; apc2=[]; apv1=[]; apv2=[]; valid=[]; n=0;
while 1
    str=fgetl(fd); if ~isstr(str), break, end
    label=SubStr(str,61,20);
    if findstr(label,'START OF ANTENNA')
        n=n+1;
        types{n}='';
        apc1(n,:)=zeros(1,3);
        apc2(n,:)=zeros(1,3);
        apv1(n,:,:)=zeros(73,19);
        apv2(n,:,:)=zeros(73,19);
        valid(n,:)=[0,inf];
        freq=0;
    elseif n<=0
        ;
    elseif findstr(label,'TYPE / SERIAL NO')
        ants{n}=strrep(strrep(SubStr(str,1,20),'NONE',''),' ',''); % omit radome-none
        if findstr(ants{n},'BLOCK')&SubStr(str,21,1)=='G'
            ants{n}=['GPS',deblank(SubStr(str,22,19))]; % gps satellite antenna
        elseif findstr(ants{n},'GLONASS')&SubStr(str,21,1)=='R'
            ants{n}=['GLO',deblank(SubStr(str,22,19))]; % glonass satellite antenna
        end
    elseif findstr(label,'VALID FROM')
        valid(n,1)=caltomjd(StrToEpoch(str,1,43)');
    elseif findstr(label,'VALID UNTIL')
        valid(n,2)=caltomjd(StrToEpoch(str,1,43)');
    elseif findstr(label,'START OF FREQUENCY')
        if strncmp(ants{n},'GLO',3)
            freq=sscanf(SubStr(str,4,3),'R%d');
        else
            freq=sscanf(SubStr(str,4,3),'G%d');
        end
    elseif findstr(label,'END OF FREQUENCY')
        freq=0;
    elseif freq==1
        if findstr(label,'NORTH / EAST / UP')
            for m=1:3, apc1(n,m)=StrToNum(str,(m-1)*10+1,10)*1E-3; end
        elseif findstr(SubStr(str,1,8),'NOAZI')
            v=LineToNum(str(9:end))*1E-3;
            for k=1:73, apv1(n,k,:)=v; end
        else
            k=round(str2num(SubStr(str,1,8))/5);
            if 0<=k&k<73, apv1(n,k+1,:)=LineToNum(str(9:end))*1E-3; end
        end
    elseif freq==2
        if findstr(label,'NORTH / EAST / UP')
            for m=1:3, apc2(n,m)=StrToNum(str,(m-1)*10+1,10)*1E-3; end
        elseif findstr(SubStr(str,1,8),'NOAZI')
            v=LineToNum(str(9:end))*1E-3;
            for k=1:73, apv2(n,k,:)=v; end
        else
            k=round(str2num(SubStr(str,1,8))/5);
            if 0<=k&k<73, apv2(n,k+1,:)=LineToNum(str(9:end))*1E-3; end
        end
    end
end

% substring --------------------------------------------------------------------
function sstr=SubStr(str,pos,ns), sstr=str(pos:min(pos+ns-1,length(str)));

% string to number -------------------------------------------------------------
function num=StrToNum(str,pos,ns)
num=sscanf(str(pos:min(pos+ns-1,length(str))),'%f',1);
if isempty(num), num=0; end

% line to number ---------------------------------------------------------------
function num=LineToNum(str)
num=zeros(1,19);
n=sscanf(str,'%f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f');
if length(n)>1, num(1:length(n))=n; num(length(n)+1:end)=n(end); end

% string to epoch --------------------------------------------------------------
function epoch=StrToEpoch(str,pos,ns)
epoch=sscanf(SubStr(str,pos,ns),'%d %d %d %d %d %f');
