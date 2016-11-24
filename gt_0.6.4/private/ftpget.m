function [stat,msg,log,bytes]=ftpget(host,login,passwd,file,ldir,redown,unzip,proto,proxy)
%-------------------------------------------------------------------------------
% [system] : GpsTools
% [module] : download files by ftp
% [func]   : download files by ftp and extract in case of compressed file
% [argin]  : host  = remote host
%            login = login user
%            passwd= password
%            file  = remote file path
%            ldir  = local directory
%           (redown)= download existing file (1:on,0:off) (default:1)
%           (unzip) = unzip/uncompact download file (1:on,0:off) (default:1)
%           (proto) = protocol ('ftp':ftp,'http':http) (default:ftp)
%           (proxy) = proxy server ([user:passwd@]server:port) (default:none)
% [argout] : stat  = result
%                 0 : ok
%                 9 : already exist
%                -1 : host/login/no directory/no file error
%                -2 : ftp exec error
%                -3 : ftp timeout
%                -4 : not supported protocol
%            msg   = result message
%            log   = ftp log
%            bytes = download file size (bytes)
% [note]   : windows
% [version]: $Revision: 20 $ $Date: 2009-05-01 04:15:33 +0900 (金, 01 5 2009) $
%            Copyright(c) 2004-2008 by T.Takasu, all rights reserved
% [history]: 04/11/19   0.1  new
%            04/12/03   0.2  fix bug on change directories
%            05/03/02   0.3  add argin redown
%            05/03/22   0.4  separate uncompact.m, add argin unzip
%            05/06/21   0.5  add status code -3 (timeout)
%                            fix bug checking existing files
%            05/07/29   0.6  use wget instead of ftp to support pasv mode
%            08/12/05   0.7  add argin proto, proxy (gt_0.6.4)
%                            support download by http
%                            support download via proxy
%                            set timeout option to 30s
%-------------------------------------------------------------------------------
if nargin<6, redown=1; end
if nargin<7, unzip=1; end
if nargin<8, proto='ftp'; end
if nargin<9, proxy=''; end
log=''; bytes=0;
if isblank(host), stat=-2; msg='Unknown host'; return, end
if isblank(login)|isblank(passwd), stat=-2; msg='Login/Passwd error'; return, end
if isblank(file), stat=-1; msg='Remote file error'; return, end
[d,f,ext]=fileparts(file);
path=fullfile(ldir,[f,ext]);
if ~redown
    if exist(path), stat=9; msg='exist'; return, end
    if any(strcmp(ext,{'.Z','.gz'}))
        [d,ff,ext]=fileparts(f);
        if exist(fullfile(ldir,f)), stat=9; msg='exist'; return, end
        if length(ext)==4&'0'<=ext(2)&ext(2)<='9'&'0'<=ext(3)&ext(3)<='9'&...
           (lower(ext(4))=='d'|lower(ext(4))=='o')
            if exist(fullfile(ldir,[ff,ext(1:3),'d'])), stat=9; msg='exist'; return, end
            if exist(fullfile(ldir,[ff,ext(1:3),'o'])), stat=9; msg='exist'; return, end
        end
    end
end
[stat,msg,bytes]=ExecFtp(host,login,passwd,ldir,file,path,proto,proxy);
if stat==0&unzip, [path,stat,msg]=uncompact(path); end
if stat==0, msg='OK'; end

% execute ftp ------------------------------------------------------------------
function [stat,msg,bytes]=ExecFtp(host,login,passwd,ldir,file,path,proto,proxy)
errs={'Unknown host','Login incorrect','Not connected','No such file',...
      'No such directory','Not a regular','File not found','Not Found',...
      'Read error','Failed to open file','Error in server'};
to=30; pre=''; bytes=0;
[dirs,f]=fileparts(which(mfilename));
cmd=fullfile(dirs,'wget');
switch proto
case 'ftp'
    if ~isempty(proxy)
       to=to*2; pre=['set ftp_proxy=http://',proxy,'/ & '];
    end
    opt=sprintf('--ftp-user=%s --ftp-password=%s --glob=off --passive-ftp -t 1 -T %d -O %s',...
                login,passwd,to,path);
    [stat,log]=dos([pre,'"',cmd,'" ',opt,' ftp://',host,file]);
case 'http'
    if ~isempty(proxy)
        to=to*2; pre=['set http_proxy=http://',proxy,'/ & '];
    end
    opt=sprintf('-t 1 -T %d -O %s',to,path);
    [stat,log]=dos([pre,'"',cmd,'" ',opt,'  http://',host,file]);
otherwise
    stat=-4; msg='Not supported protocol'; return
end
for n=1:length(errs)
   if findstr(lower(log),lower(errs{n}))
       if exist(path), delete(path), end
       msg=errs{n}; stat=-1; return
   end
end
if stat~=0
    if ~exist(path), stat=-2; msg='Ftp/http exec error';
    else stat=-3; msg='Ftp/http get error'; delete(path); end
else
    info=dir(path); bytes=info.bytes;
    if bytes>0, msg='OK'; else stat=-3; msg='Empty file'; delete(path); end
end

% check blank/number string ----------------------------------------------------
function stat=isblank(str)
stat=min(isspace(str)); if isempty(stat), stat=1; end

function stat=isnumber(str), stat=all('0'<=str&str<='9');
