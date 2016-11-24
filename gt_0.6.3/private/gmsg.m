function abort=gmsg(varargin)
%-------------------------------------------------------------------------------
% [system] : GpsTools
% [module] : display message
% [func]   : display message and show progress
% [argin]  : options
%               'seth',h  : set output handle
%                prog     : progress (0-1)
%                msg      : display messages
% [argout] : abort = abort button pusshed
% [note]   :
% [version]: $Revision: 5 $ $Date: 06/07/14 16:12 $
%            Copyright(c) 2004-2006 by T.Takasu, all rights reserved
% [history]: 04/12/06  0.1  new
%            05/03/19  0.2  add function setting output handles
%-------------------------------------------------------------------------------
persistent hs ts prog, if isempty(prog), prog=0; end
abort=0;
if ischar(varargin{1})
    if strcmp(varargin{1},'seth')
        hs=varargin{2};
    else
        msg=sprintf(varargin{:});
        if isempty(hs), disp(msg); else set(hs(1),'string',msg); end
    end
elseif ~isempty(hs)
    prog=varargin{1}; gut('setprogh',hs(2),prog);
end
if ~isempty(hs), abort=get(hs(3),'userdata'); end

if ticks, h=gcf; drawnow; figure(h); end

% check tick -------------------------------------------------------------------
function f=ticks
persistent tp, ts=now; if 86400*(ts-tp)<0.5, f=0; else f=1; tp=ts; end
