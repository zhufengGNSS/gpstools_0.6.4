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
% [version]: $Revision: 3 $ $Date: 06/07/08 1:16 $
%            Copyright(c) 2004-2006 by T.Takasu, all rights reserved
% [history]: 06/03/26  0.1  separated from gpsestd.m
%            06/06/24  0.2  modify specifications
%-------------------------------------------------------------------------------
persistent flg logs, if isempty(logs), logs={}; end
log={};
switch varargin{1}
case 'start', flg=1; logs={};
case 'stop',  flg=0; log=logs; logs={};
otherwise
    msg=sprintf(varargin{:}); if flg, logs={logs{:},msg}; else disp(msg); end
end
