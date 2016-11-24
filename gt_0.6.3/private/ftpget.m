function [stat,msg,log,bytes]=ftpget(host,login,passwd,file,ldir,redown,unzip)
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
% [argout] : stat  = result
%                 1 : alreay exist
%                 0 : ok
%                -1 : host/login/no directory/no file error
%                -2 : ftp exec error
%                -3 : ftp timeout
%            msg   = result message
%            log   = ftp log
%            bytes = download file size (bytes)
% [note]   : windows
% [version]: $Revision: 2 $ $Date: 06/07/08 1:16 $
%            Copyright(c) 2004-2006 by T.Takasu, all rights reserved
% [history]: 04/11/19   0.1  new
%            04/12/03   0.2  fix bug on change directories
%            05/03/02   0.3  add argin redown
%            05/03/22   0.4  separate uncompact.m, add argin unzip
%            05/06/21   0.5  add status code -3 (timeout)
%                            fix bug checking existing files
%-------------------------------------------------------------------------------
if nargin<6, redown=1; end
if nargin<7, unzip=1; end
log=''; bytes=0;
if isblank(host), stat=-2; msg='Unknown host'; return, end
if isblank(login)|isblank(passwd), stat=-2; msg='Login/Passwd error'; return, end
if isblank(file), stat=-1; msg='Remote file error'; return, end
[d,f,ext]=fileparts(file);
path=fullfile(ldir,[f,ext]);
if ~redown
    if exist(path), stat=1; msg='exist'; return, end
    if any(strcmp(ext,{'.Z','.gz'}))
        [d,ff,ext]=fileparts(f);
        if exist(fullfile(ldir,f)), stat=1; msg='exist'; return, end
        if length(ext)==4&'0'<=ext(2)&ext(2)<='9'&'0'<=ext(3)&ext(3)<='9'&...
           (lower(ext(4))=='d'|lower(ext(4))=='o')
            if exist(fullfile(ldir,[ff,ext(1:3),'d'])), stat=1; msg='exist'; return, end
            if exist(fullfile(ldir,[ff,ext(1:3),'o'])), stat=1; msg='exist'; return, end
        end
    end
end
[stat,msg,bytes]=ExecFtp(host,login,passwd,ldir,file,path);
if stat==0&unzip, [path,stat,msg]=uncompact(path); end
if stat==0, msg='OK'; end

% execute ftp ------------------------------------------------------------------
function [stat,msg,bytes]=ExecFtp(host,login,passwd,ldir,file,path)
errs={'Unknown host','Login incorrect','Not connected','No such file',...
      'Not a regular','File not found','Failed to open file'};
bytes=0;
fid=fopen('ftp.tmp','wt');
if fid==-1, stat=-2; msg='Ftp ctr error'; return, end
fprintf(fid,'user %s %s\nbin\nlcd %s\nget %s\nquit\n',login,passwd,ldir,file);
fclose(fid);
[stat,log]=dos(['ftp -v -n -i -s:ftp.tmp ',host]);
delete('ftp.tmp')
if stat~=0, msg='Ftp execution error'; return, end
for n=1:length(errs)
   if findstr(log,errs{n})
       if exist(path), delete(path), end
       msg=errs{n}; stat=-1; return
   end
end
info=dir(path);
if isempty(info), stat=-3; msg='Ftp get error'; return, end
bytes=info.bytes;
if bytes==0, delete(path), stat=-3; msg='Empty file'; return, end
stat=0; msg='OK';

% check blank/number string ----------------------------------------------------
function stat=isblank(str)
stat=min(isspace(str)); if isempty(stat), stat=1; end

function stat=isnumber(str), stat=all('0'<=str&str<='9');
