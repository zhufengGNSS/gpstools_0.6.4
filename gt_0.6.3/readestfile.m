function [td,time,xs,vs,prm,type,fb,name]=readestfile(files)
%-------------------------------------------------------------------------------
% [system] : GpsTools
% [module] : read estimation result files
% [func]   : read estimation result files
% [argin]  : files = estimation result files
% [argout] : td    = first date (mjd)
%            time  = time vector (sec)
%            xs    = estimation results
%            vs    = estimation variences
%            prm   = estimation parameters (first)
%            type  = type
%            fb    = forwart/back
%            name  = satellite/station name
% [note]   :
% [version]: $Revision: 1 $ $Date: 06/07/08 14:21 $
%            Copyright(c) 2004-2006 by T.Takasu, all rights reserved
% [history]: 04/11/17  0.1  new
%            04/12/08  0.2  argin tunit in hr instead of sec
%-------------------------------------------------------------------------------
td=[]; time=[]; xs=[]; vs=[]; prm=[]; type=''; fb=''; name='';

if ~iscell(files), files={files}; end

for n=1:length(files)
    [path,file,ext]=fileparts(files{n});
    if strcmp(ext,'.mat')
        [typ,file]=strtok(file,'_');
        [nam,file]=strtok(file,'_');
        if (isempty(type)|strcmp(typ,type))&(isempty(name)|strcmp(nam,name))
            type=typ;
            name=nam;
            [epoch,t,data,covs,pr]=loadestfile(files{n});
            if ~isempty(epoch)
                [tdd,ts]=caltomjd(epoch);
                if isempty(td), td=tdd; prm=pr; end
                time=[time;t+ts+(tdd-td)*86400];
                xs=[xs;data];
                if ~isempty(covs), vs=[vs;covs]; else vs=[vs;repmat(nan,size(xs))]; end
            end
        else
            disp(['can not marged : ',files{n}])
        end
    end
end
if ~isempty(time)
    [time,i]=sortrows(time);
    xs=xs(i,:);
    vs=vs(i,:);
    dd=floor(time(1)/86400);
    td=td+dd;
    time=time-dd*86400;
end
fb=type(end);
type=type(1:end-1);

% load estimate file -------------------------------------------------------------------------------
function [epoch,time,data,covs,prm]=loadestfile(file)
epoch=[]; time=[]; data=[]; covs=[]; prm=[];
if exist(file), load(file), end
