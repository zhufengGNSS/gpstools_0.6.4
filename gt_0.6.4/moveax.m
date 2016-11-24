function moveax(ax)
% Add move/zoom function to axes. To move x-axis, push mouse left-button at
% lower-edge of axis and drag. To move y-axis, push it at left-edge and drag.
% To zoom/unzoom axes, use mouse right-button. moveax effects to current axes.
% To specify axes handle h, use moveax(h).
%
% ver.1.0 2006/10/13 by T.T

if nargin<1, ax=gca; elseif strncmp(ax,'cb_',3), feval(ax); return; end

p.mov=0; p.pnt=[0,0]; p.xl=[0,1]; p.yl=[0,1];
set(ax,'userdata',p);
set(get(ax,'parent'),'windowbuttondownfcn',[mfilename,' cb_btndown'],...
    'windowbuttonmotionfcn',[mfilename,' cb_btndrag'],...
    'windowbuttonupfcn',[mfilename,' cb_btnup']);

function cb_btndown
ax=get(gcf,'currentaxes'); if isempty(ax), return; end
p=get(ax,'userdata'); if ~isstruct(p)|~isfield(p,'mov'), return; end
siz=get(gcf,'position');
pos=get(ax,'position');
pos=[siz(3:4).*pos(1:2),siz(3:4).*(pos(1:2)+pos(3:4))];
p.pnt=get(gcf,'currentpoint');
if     abs(p.pnt(2)-pos(2))<10&pos(1)<p.pnt(1)&p.pnt(1)<pos(3), p.mov=1;
elseif abs(p.pnt(1)-pos(1))<10&pos(2)<p.pnt(2)&p.pnt(2)<pos(4), p.mov=2;
else return; end
if strcmp(get(gcf,'selectiontype'),'alt'), p.mov=p.mov+2; end % right button
p.xl=xlim; p.yl=ylim;
set(ax,'userdata',p);
cur={'right','top'};
if p.mov<=2, set(gcf,'pointer',cur{p.mov}); else setcurzoom(p.mov==3); end

function cb_btndrag
ax=get(gcf,'currentaxes'); if isempty(ax), return; end
p=get(ax,'userdata'); if ~isstruct(p)|~isfield(p,'mov')|p.mov<=0, return; end
siz=get(gcf,'position');
pos=get(ax,'position');
dp=(get(gcf,'currentpoint')-p.pnt)./siz(3:4)./pos(3:4);
switch p.mov
case 1, p.xl=p.xl-(p.xl(2)-p.xl(1))*dp(1); % move x axis
case 2, p.yl=p.yl-(p.yl(2)-p.yl(1))*dp(2); % move y axis
case 3, p.xl=p.xl(1)+(p.xl(2)-p.xl(1))*(0.5+[-0.5,0.5]*0.1^dp(1)); % zoom x axis
case 4, p.yl=p.yl(1)+(p.yl(2)-p.yl(1))*(0.5+[-0.5,0.5]*0.1^dp(2)); % zoom y axis
end
xlim(p.xl); ylim(p.yl);

function cb_btnup
ax=get(gcf,'currentaxes'); if isempty(ax), return; end
p=get(ax,'userdata'); if ~isstruct(p)|~isfield(p,'mov')|p.mov<=0, return; end
p.mov=0;
set(ax,'userdata',p);
set(gcf,'pointer','arrow');

function setcurzoom(opt)
o=nan;
p=[o,o,o,o,o,o,o,2,1,2,o,o,o,o,o,o
   o,o,o,o,o,o,2,1,1,1,2,o,o,o,o,o
   o,o,o,o,o,2,1,1,1,1,1,2,o,o,o,o
   o,o,o,o,2,1,1,1,1,1,1,1,2,o,o,o
   o,o,o,2,1,1,1,1,1,1,1,1,1,2,o,o
   o,o,o,2,2,2,2,2,1,2,2,2,2,2,o,o
   o,o,o,o,o,o,o,2,1,2,o,o,o,o,o,o
   o,o,o,o,o,o,o,2,1,2,o,o,o,o,o,o
   o,o,o,o,o,2,2,2,1,2,2,2,o,o,o,o
   o,o,o,o,o,2,1,1,1,1,1,2,o,o,o,o
   o,o,o,o,o,2,2,2,1,2,2,2,o,o,o,o
   o,o,o,o,o,o,o,2,1,2,o,o,o,o,o,o
   o,o,o,o,o,2,2,2,1,2,2,2,o,o,o,o
   o,o,o,o,o,2,1,1,1,1,1,2,o,o,o,o
   o,o,o,o,o,o,2,1,1,1,2,o,o,o,o,o
   o,o,o,o,o,o,o,2,1,2,o,o,o,o,o,o];
if opt, p=fliplr(p'); end
set(gcf,'pointer','custom','pointershapecdata',p,'pointershapehotspot',[9,9])
