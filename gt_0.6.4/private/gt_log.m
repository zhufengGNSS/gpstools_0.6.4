function log=gt_log(varargin)
%-------------------------------------------------------------------------------
% [system] : GpsTools
% [module] : write log
% [func]   : write log
% [argin]  : msg,... = log message
%                      msg{1}     : command or sprintf message format
%                          msg{1}='start' : start recording log
%                          msg{1}='stop'  : stop recording log and get log
%                      msg{2:end} : sprintf data
% [argout] : log     = recorded log messages (cell array)
% [note]   : if not recording mode log output to console
% [version]: $Revision: 16 $ $Date: 2008-12-12 15:49:30 +0900 (金, 12 12 2008) $
%            Copyright(c) 2004-2006 by T.Takasu, all rights reserved
% [history]: 06/03/26  0.1  separated from gpsestd.m
%            06/06/24  0.2  modify specifications
%            08/11/26  0.3  suppress duplicated messages (gt_0.6.4)
%-------------------------------------------------------------------------------
persistent flg logs, if isempty(logs), logs={}; end
log={};
switch varargin{1}
case 'start', flg=1; logs={};
case 'stop',  flg=0; log=logs; logs={};
otherwise
    msg=sprintf(varargin{:});
    if flg
        if isempty(logs)|~strcmp(msg,logs{end}), logs={logs{:},msg}; end
    else
        disp(msg);
    end
end
