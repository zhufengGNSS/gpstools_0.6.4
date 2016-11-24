function prm=loadprm(prmname,prmdefault)
%-------------------------------------------------------------------------------
% [system] : GpsTools
% [module] : load parameters
% [func]   : load parameters from the following precedence in search pathes
%               (1) from <prmname> (mat-file) in search path
%               (2) from <prmname> (m-file) in search path
%               (3) from <gtpath>/settings/<prmname>.mat
%               (4) from <gtpath>/settings/<prmdefault>.m
%               (<gtpath> : GT commnad path directory)
% [argin]  : prmname    = parameter
%           (prmdefault)= default parameter
% [argout] : prm = parameter struct ([]:no parameter)
% [note]   : parameter mat-file must contain variable 'prm'
%            parameter m-file must return parameter struct
% [version]: $Revision: 12 $ $Date: 2008-11-25 10:02:15 +0900 (火, 25 11 2008) $
%            Copyright(c) 2004-2006 by T.Takasu, all rights reserved
% [history]: 04/12/11  0.1  new
%            06/02/13  0.2  set hierarchical fields to parameter struct
%            06/06/30  0.3  modify parameter file position
%-------------------------------------------------------------------------------
if nargin<2, prmdefault=''; end
prm=[]; fields={};

% default settings directory path
[d,f]=fileparts(which(mfilename)); sdir=[d,filesep,'settings'];

if exist(prmname)==2
    [dirs,file,ext]=fileparts(prmname);
    if isempty(ext)
        [dirs,file,ext]=fileparts(which(prmname));
    end
    if strcmp(ext,'.mat')
        load(prmname,'prm');
    elseif strcmp(ext,'.m')
        wd=cd; cd(dirs); prm=feval(file); cd(wd);
    end
else
    prmfile=fullfile(sdir,[prmname,'.mat']);
    if exist(prmfile)==2
        load(prmfile,'prm');
    end
end

% use default parameter
if ~isempty(prmdefault)
    wd=cd; cd(sdir); prm=copyfield(feval(prmdefault),prm); cd(wd);
end

% copy fields ------------------------------------------------------------------
function prms=copyfield(prms,prmd)
if isstruct(prms)
    fnames=fieldnames(prms);
    for n=1:length(fnames)
        if isfield(prmd,fnames{n})
            field=copyfield(getfield(prms,fnames{n}),getfield(prmd,fnames{n}));
            prms=setfield(prms,fnames{n},field);
        end
    end
elseif iscell(prms)|ischar(prms)|size(prms,1)<size(prmd,1)|size(prms,2)<size(prmd,2)
    prms=prmd;
else
    m=size(prmd,1); n=size(prmd,2); prms(1:m,1:n)=prmd(1:m,1:n);
end
