function [disps,phass]=readoload(rcvs,file)
%-------------------------------------------------------------------------------
% [system] : GpsTools
% [module] : read ocean loading parameters
% [func]   : read ocean loading parameters
% [argin]  : rcvs  = station list
%            file  = ocean loading parameters file path
% [argout] : disps = ocean loading amplitudes (m) [radial,east,north]
%                  disps(n,:,m) = const.n,rcvs{m} amplitude
%            phass = ocean loading phase (deg) [radial,east,north]
%                  covs(n,:,m) = const.n,rcvs{n} phase
% [note]   : constituents : M2,S2,N2,K2,K1,O1,P1,Q1,MF,MM,SSA
% [version]: $Revision: 16 $ $Date: 2008-12-12 15:49:30 +0900 (金, 12 12 2008) $
%            Copyright(c) 2004-2006 by T.Takasu, all rights reserved
% [history]: 04/11/02  0.1  new
%            04/12/08  0.2  change argin oloadfile
%            05/01/05  0.3  support case of no ocean loading files
%            06/06/25  0.4  add warning messages
%            08/11/26  0.5  fix bug on error-stop if format error
%-------------------------------------------------------------------------------
if ischar(rcvs), rcvs={rcvs}; end, rcvs=upper(rcvs);
disps=repmat(nan,[11,3,length(rcvs)]); phass=disps;
if isempty(file), return, end

[rrcvs,rdisps,rphass]=ReadOloadFile(file);
if isempty(rrcvs), return; end

for n=1:length(rcvs)
    i=find(strcmp(rcvs{n},rrcvs));
    if ~isempty(i)
        disps(:,:,n)=rdisps(:,:,i); phass(:,:,n)=rphass(:,:,i);
    else
        gt_log('no ocean loading prms   : %s file=%s',rcvs{n},file);
    end
end

% read blq ocean loading parameters file ---------------------------------------
function [rcvs,disps,phass]=ReadOloadFile(file)
rcvs={}; disps=[]; phass=[];
fd=fopen(file,'rt');
if fd<0, gt_log('no ocean loading file   : %s',file); return; end
[rcvs,disps,phass]=ReadBlqFile(fd);
fclose(fd);

% read blq ocean loading parameters --------------------------------------------
function [rcvs,disps,phass]=ReadBlqFile(fd)
start=0; n=0; m=0; rcvs={}; disps=[]; phass=[];
while 1
    str=fgetl(fd); if ~isstr(str), break, end
    if strfind(str,'$$ END TABLE'), break
    elseif strcmp(SubStr(str,1,2),'$$') ; % comment
    elseif m<=0|7<=m
        rcvs={rcvs{:},deblank(str(3:end))}; n=n+1; m=1;
    elseif m<4
        for k=1:11, disps(k,m,n)=StrToNum(str,7*k-4,6); end, m=m+1;
    elseif m<7
        for k=1:11, phass(k,m-3,n)=StrToNum(str,7*k-4,6); end, m=m+1;
    end
end

% substring --------------------------------------------------------------------
function sstr=SubStr(str,pos,ns), sstr=str(pos:min(pos+ns-1,length(str)));

% string->number ---------------------------------------------------------------
function num=StrToNum(str,pos,ns)
num=sscanf(str(pos:min(pos+ns-1,length(str))),'%f',1);
if isempty(num), num=nan; end
