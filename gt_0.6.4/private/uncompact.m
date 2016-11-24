function [path,stat,msg]=uncompact(path,varargin)
%-------------------------------------------------------------------------------
% [system] : GpsTools
% [module] : unzip/uncompact file
% [func]   : unzip and uncompact file
% [argin]  : path  = file path
%            (opt) = options
%                'dir',dirs : output directory (default:same as original file)
%                'org'      : retain original file (default:delete orginal file)
% [argout] : path  = unzip/uncompact file path
%            stat  = status (0:ok,other:error)
%            msg   = error message
% [note]   :
% [version]: $Revision: 12 $ $Date: 2008-11-25 10:02:15 +0900 (火, 25 11 2008) $
%            Copyright(c) 2004-2006 by T.Takasu, all rights reserved
% [history]: 05/03/19   0.1  separated from ftpget
%            05/04/25   0.2  add argin opt
%            06/07/02   0.3  fix bug of crx2rnx return code
%-------------------------------------------------------------------------------
stat=0; msg='OK'; n=1; odir=''; org=0;
while n<=nargin-1
    switch varargin{n}
    case 'dir', odir=varargin{n+1}; n=n+2;
    case 'org', org=1; n=n+1;
    otherwise, error('argin error'), end
end
[dirs,file,ext]=fileparts(path);
gzip=which('gzip.exe');
if ~isempty(gzip)&(strcmp(ext,'.Z')|strcmp(ext,'.gz'))
    if ~isempty(odir)|org
        if ~isempty(odir), dirs=odir; end
        [stat,msg]=dos(['"',gzip,'" -f -d -c "',path,'"',' >"',fullfile(dirs,file),'"']);
        if stat~=0, delete(fullfile(dirs,file)); return, end
        org=0;
    else
        [stat,msg]=dos(['"',gzip,'" -f -d "',path,'"']);
    end
    if stat~=0, return, end
    path=fullfile(dirs,file);
end
[dirs,file,ext]=fileparts(path);
crx2rnx=which('crx2rnx.exe');
if ~isempty(crx2rnx)&length(ext)==4&...
   all('0'<=ext(2:3)&ext(2:3)<='9')&lower(ext(4))=='d'
    if ext(4)=='d', ext(4)='o'; else ext(4)='O'; end
    npath=fullfile(dirs,[file,ext]);
    [stat,msg]=dos(['"',crx2rnx,'" < "',path,'" > "',npath,'"']);
    if stat~=0
        if exist(npath), delete(npath); end
        return
    end
    if ~org, delete(path); end
    path=npath;
end
