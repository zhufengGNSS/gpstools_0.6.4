function h=t_uieditbug(varargin)
% test uicontrol bug of matlab

varargin{:}

if nargin>=1&strncmp(varargin{1},'cb',2)
    feval(varargin{:}); return;
end

set(0,'defaultuicontrolfontname','times new roman');

figure

h=uicontrol('style','edit','position',[10 10 100 20]);

set(gcf,'keypressfcn',[mfilename,' cbkeyFig']);
set(h,'keypressfcn',[mfilename,' cbkeyEdit']);

% callback
function cbkeyEdit
disp('key pressed in edit');

function cbkeyFig
disp('key pressed in figure');
