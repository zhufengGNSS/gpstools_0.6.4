function path=gfilepath(path,file,day,name,opt)
%-------------------------------------------------------------------------------
% [system] : GpsTools
% [module] : extend file path
% [func]   : replace keywords and extend file path
% [argin]  : path  = directory path
%           (file) = file name (default:'')
%           (day)  = day [year,month,day] (default:[])
%           (name) = satellite/station name (default:'')
%           (opt)  = option
%                    0 : not make directory (default)
%                    1 : make directory if not exist
% [argout] : path  = file path
% [note]   : keywords in file path are replaced as follows
%               %S -> satellite/station name (UPPER case)
%               %s -> satellite/station name (lower case)
%               %G -> GSI station code (tail 4 characters)
%               %Y -> year           (1990-2010)
%               %y -> year           (00-99)
%               %m -> month          (01-12)
%               %d -> day            (01-31)
%               %n -> day of year    (001-366)
%               %W -> gpsweek no     (0000-9999)
%               %D -> day of gpsweek (0-6)
%               %M -> iers bulletin month (1- )
%               %P -> install path
%               %RO-> rinex obs      (%s%n0.%yo)
%               %RD-> rinex cobs     (%s%n0.%yd)
%               %RN-> rinex nav (GPS)(%s%n0.%yn)
%               %RG-> rinex nav (GLONASS) (%s%n0.%yg)
%               %RM-> rinex met      (%s%n0.%ym)
%               %RS-> rinex sum      (%s%n0.%ys)
% [version]: $Revision: 12 $ $Date: 2008-11-25 10:02:15 +0900 (火, 25 11 2008) $
%            Copyright(c) 2004-2006 by T.Takasu, all rights reserved
% [history]: 05/02/28  0.1  separated from gpsestd.m
%            05/05/15  0.2  fix bug that dirs does not made if file empty
%            06/07/01  0.3  improve performance, delete replacement of %EA
%-------------------------------------------------------------------------------
if nargin<2, file=''; end
if nargin<3, day=[]; end
if nargin<4, name=''; end
if nargin<5, opt=0; end
s=filesep;
if isempty(path)|path(end)==s, path=[path,file]; else path=[path,s,file]; end
if findstr(path,'%')
    if findstr(path,'%P')
        [dirs,f]=fileparts(which(mfilename)); [dirs,f]=fileparts(dirs);
        path=strrep(path,'%P',dirs);
    end
    if findstr(path,'%R')
        path=strrep(path,'%RO','%s%n0.%yo');
        path=strrep(path,'%RN','%s%n0.%yn');
        path=strrep(path,'%RG','%s%n0.%yg');
        path=strrep(path,'%RM','%s%n0.%ym');
        path=strrep(path,'%RD','%s%n0.%yd');
        path=strrep(path,'%RS','%s%n0.%ys');
    end
    if ~isempty(name)
        name=char(name);
        if findstr(path,'%s'), path=strrep(path,'%s',lower(name)); end
        if findstr(path,'%S'), path=strrep(path,'%S',upper(name)); end
        if findstr(path,'%G'), path=strrep(path,'%G',name(max(1,end-3):end)); end
    end
    if ~isempty(day)
        dn=datenum(day); d0=723186; % datenum(1980,1,6);
        if findstr(path,'%Y'), path=strrep(path,'%Y',sprintf('%04d',day(1))); end
        if findstr(path,'%y'), path=strrep(path,'%y',sprintf('%02d',mod(day(1),100))); end
        if findstr(path,'%m'), path=strrep(path,'%m',sprintf('%02d',day(2))); end
        if findstr(path,'%d'), path=strrep(path,'%d',sprintf('%02d',day(3))); end
        if findstr(path,'%n'), path=strrep(path,'%n',sprintf('%03d',floor(dn-datenum(day(1),1,1))+1)); end
        if findstr(path,'%W'), path=strrep(path,'%W',sprintf('%03d',floor((dn-d0)/7))); end
        if findstr(path,'%D'), path=strrep(path,'%D',sprintf('%d',floor(mod(dn-d0,7)))); end
        if findstr(path,'%M'), path=strrep(path,'%M',sprintf('%d',(day(1)-1988)*12+day(2)-1)); end
    end
end
if opt
    if isempty(file) dirs=path; else [dirs,file]=fileparts(path); end
    [root,dirs]=strtok(dirs,filesep);
    if ~isempty([root,dirs])&~exist([root,dirs])
        mkdir(root,dirs);
    end
end
