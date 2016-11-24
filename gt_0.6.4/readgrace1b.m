function [epoch,time,index,data,anc,type,info]=readgrace1b(file)
%-------------------------------------------------------------------------------
% [system] : GpsTools
% [module] : read GRACE Level 1B data file
% [func]   : read GRACE Level 1B data file
% [argin]  : file  = file path
% [argout] : epoch = epoch [y,m,d,h,m,s]
%            time  = time (sec)
%            index = data index
%            data  = data
%            anc   = ancylary data
%            type  = data type
%            info  = header info
% [note]   :
% [version]: $Revision: 12 $ $Date: 2008-11-25 10:02:15 +0900 (火, 25 11 2008) $
%            Copyright(c) 2004-2006 by T.Takasu, all rights reserved
% [history]: 05/06/09   0.1  new
%-------------------------------------------------------------------------------
epoch=[]; time=[]; index=[]; data=[]; anc=[]; type=''; info={};
fd=fopen(file,'r','ieee-be');
if fd<0, disp(['warning : grace level 1b file open error : ',file]), return, end
[p,f,e]=fileparts(file); type=f(1:5);
epoch=[str2num(f(7:10)),str2num(f(12:13)),str2num(f(15:16)),0,0,0];
t0=(caltomjd(epoch)-caltomjd([2000,1,1,12,0,0,0]))*86400;
info=ReadHeader(fd);
switch type
case 'CLK1B', [time,index,data,anc]=ReadSatClk(fd,t0);
case 'GPS1B', [time,index,data,anc]=ReadGpsObs(fd,t0);
case 'GNV1B', [time,index,data,anc]=ReadSatEph(fd,t0);
case {'VGN1B','VGO1B','VGB1B','VCM1B','VKB1B','VSL1b'}
    [time,index,data,anc]=ReadGpsApc(fd,t0);
otherwise
    disp(['warning : no type supported : ',type])
end
fclose(fd);

% read header ------------------------------------------------------------------
function info=ReadHeader(fd)
info={};
while 1
    str=fgetl(fd); if ~isstr(str), break, end
    info={info{:},deblank(str)};
    if ~isempty(findstr(str,'END OF HEADER')), break, end
end

% read satellite clock ---------------------------------------------------------
function [time,index,data,sig]=ReadSatClk(fd,t0)
time=[]; index=[]; data=[]; sig=[]; n=0;
while 1
    t=fread(fd,1,'int32')-t0;         % rcv_time
    s=fread(fd,2,'uchar');            % GRACE_id/clock_id
    d=fread(fd,4,'double');           % eps_time/eps_err/eps_drift/drift_err
    q=fread(fd,1,'uchar');            % qualflg
    if length(d)<4, break
    elseif q==0
       n=n+1;
       time(n,1)=t;
       index(n,:)=[s(1)-'A'+1,s(2)];
       data(n,:)=d([1,3])';
       sig(n,:)=d([2,4])';
    end
end

% read gps observation data ----------------------------------------------------
function [time,index,data,anc]=ReadGpsObs(fd,t0)
time=[]; index=[]; data=[]; anc=[]; n=0;
while 1
    t=fread(fd,1,'int32')-t0;         % rcvtime_intg
    t=t+fread(fd,1,'int32')*1E-6;     % rcvtime_frac
    s=fread(fd,3,'uchar');            % GRACE_id/prn_id/ant_id
    p=fread(fd,1,'ushort');           % prod_flag
    q=fread(fd,1,'uchar');            % qualflg
    d=fread(fd,6,'double');           % CA/L1/L2 pseudo range(m)/phase(m)
    a=fread(fd,6,'ushort');           % CA/L1/L2 SNR/CH
    if length(a)<6, break
    else
       n=n+1; disp(sprintf('n=%d,t=%.0f',n,t))
       time(n,1)=t;
       index(n,:)=[s(1)-'A'+1,s(2:3)'];
       data(n,:)=d';
       anc(n,:)=a';
    end
end

% read gps orbit ---------------------------------------------------------------
function [time,index,data,anc]=ReadSatEph(fd,t0)
time=[]; index=[]; data=[]; anc=[]; n=0;
while 1
    t=fread(fd,1,'int32')-t0;         % gps_time
    s=fread(fd,2,'uchar');            % GRACE_id/coord_ref
    d=fread(fd,3,'double');           % xpos/ypos/zpos
    a=fread(fd,3,'double');           % xpos_err/ypos_err/zpos_err
    d=[d;fread(fd,3,'double')];       % xvel/yvel/zvel
    a=[a;fread(fd,3,'double')];       % xvel_err/yvel_err/zvel_err
    q=fread(fd,1,'uchar');            % qualflg
    if length(a)<6, break
    elseif q==0
       n=n+1;
       time(n,1)=t;
       index(n,:)=[s(1)-'A'+1,s(2)];
       data(n,:)=d';
       anc(n,:)=a';
    end
end

% read antenna phase center ----------------------------------------------------
function [time,index,data,anc]=ReadGpsApc(fd,t0)
time=[]; index=[]; data=[]; anc=[]; n=0;
while 1
    t=fread(fd,1,'int32')-t0;         % gps_time
    s=fread(fd,1,'uchar');            % GRACE_id
    d=fread(fd,4,'double');           % mag/cosx/cosy/cosz
    q=fread(fd,1,'uchar');            % qualflg
    if length(d)<4, break
    else
       n=n+1;
       time(n,1)=t;
       index(n,:)=[s(1)-'A'+1,q];
       data(n,:)=d';
    end
end
