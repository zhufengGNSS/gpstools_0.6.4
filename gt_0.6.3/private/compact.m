function [path,stat,msg]=compact(path,varargin)
%-------------------------------------------------------------------------------
% [system] : GpsTools
% [module] : zip/compact file
% [func]   : zip and compact file
% [argin]  : path  = file path
%            (opt) = options
%                'dir',dir : output directory (default:same as original file)
%                'org' : keep original file (default:delete orginal file)
%                'gz'  : .gz as sufix (default:.Z)
% [argout] : path  = zip/compact file path
%            stat  = status (0:ok,other:error)
%            msg   = error message
% [note]   :
% [version]: $Revision: 3 $ $Date: 06/07/13 8:03 $
%            Copyright(c) 2004-2006 by T.Takasu, all rights reserved
% [history]: 05/04/25   0.1  new
%-------------------------------------------------------------------------------
stat=0; msg='OK'; n=1; odir=''; org=0; suff='.Z';
while n<=nargin-1
    switch varargin{n}
    case 'dir', odir=varargin{n+1}; n=n+2;
    case 'org', org=1; n=n+1;
    case 'gz',  suff='.gz'; n=n+1;
    otherwise, n=n+1; end
end
[dirs,file,ext]=fileparts(path);
rnx2crx=which('rnx2crx.exe');
if ~isempty(rnx2crx)&length(ext)==4&...
   all('0'<=ext(2:3)&ext(2:3)<='9')&lower(ext(4))=='o'
    if ext(4)=='o', ext(4)='d'; else ext(4)='D'; end
    npath=fullfile(dirs,[file,ext]);
    [stat,msg]=dos(['"',rnx2crx,'" < "',path,'" > "',npath,'"']);
    if stat~=0
        if exist(npath), delete(npath); end
        return
    end
    if ~org, delete(path); end
    path=npath;
end
[dirs,file,ext]=fileparts(path);
gzip=which('gzip.exe');
if ~isempty(gzip)
    if ~isempty(odir)|org
        if ~isempty(odir), dirs=odir; end
        npath=fullfile(dirs,[file,ext,suff]);
        [stat,msg]=dos(['"',gzip,'" -f -c "',path,'"',' >"',npath,'"']);
        if stat~=0, delete(npath); return, end
        if ~org, delete(path); end
    else
        [stat,msg]=dos(['"',gzip,'" -f -S ',suff,' "',path,'"']);
        npath=[path,suff];
    end
    if stat~=0, return, end
    path=npath;
end
