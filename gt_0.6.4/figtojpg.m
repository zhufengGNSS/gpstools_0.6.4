function figtojpg(figs,dirs)
%-------------------------------------------------------------------------------
% [system] : GpsTools
% [module] : save current figures to jpeg files
% [func]   : save current figures to jpeg files
% [argin]  : (figs) = figure numbers (default:all)
%            (dirs) = directory      (default:current)
% [argout] : none
% [note]   : Figure No.1->fig01.jpg
%            Figure No.2->fig02.jpg
%            ...
% [version]: $Revision: 12 $ $Date: 2008-11-25 10:02:15 +0900 (火, 25 11 2008) $
%            Copyright(c) 2004-2006 by T.Takasu, all rights reserved
% [history]: 04/10/03  0.1  new
%-------------------------------------------------------------------------------
if nargin<1, figs=[]; end
if nargin<2, dirs=''; end
for h=get(0,'Children')'
    if isempty(figs)|any(h==figs)
        file=fullfile(dirs,sprintf('fig%02d.jpg',h));
        c=get(h,'color'); set(h,'color','w');
        imwrite(frame2im(getframe(h)),file,'jpeg');
        set(h,'color',c);
    end
end
