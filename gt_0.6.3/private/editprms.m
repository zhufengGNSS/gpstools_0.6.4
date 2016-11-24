function [list,prms,sel]=editprms(title,list,prms,name,prm)
%-------------------------------------------------------------------------------
% [system] : GpsTools
% [module] : parameter list editor dialog
% [func]   : parameter list editor dialog
% [argin]  : title = dialog title
%            list  = parameters name list
%            prms  = parameters (cell array)
%            name  = added parameter name
%            prm   = added parameter
% [argout] : list  = parameters name list
%            prms  = parameters
%            sel   = select parameter index ([]:cancel)
% [note]   : modal dialog
% [version]: $Revision: 3 $ $Date: 06/07/08 9:36 $
%            Copyright(c) 2004-2006 by T.Takasu, all rights reserved
% [history]: 05/06/07  0.1  new
%-------------------------------------------------------------------------------
if strncmp(title,'cb',2), feval(title)
else [list,prms,sel]=EditPrms(title,list,prms,name,prm); end

% edit parameters list body ----------------------------------------------------
function [list,prms,sel]=EditPrms(title,list,prms,name,prm)
h=gut('newdlg','pdlg',title,[285,286]);
gut('newlist','list',[3,50,280,232],list,prms,2,[mfilename,' cbPrmClk']);
set(gut('newedit','prm',[4,29,278,20],name,1),'userdata',prm);
gut('newbtnh','btns',[3,3,140,23],{'Add/Set','Delete'},...
    {[mfilename,' cbSet'],[mfilename,' cbDel']});
gut('newokcancelbtn','',[143,3,140,23]);
if isempty(prm), gut('setdis',{'prm','btns_1'});
else gut('setstring','prm','newlist'); end
if ~gut('waitok'); sel=[]; return, end
[list,prms]=gut('getlist','list');
sel=min(gut('getval','list'));
close

% callback on set --------------------------------------------------------------
function cbSet
name=gut('getstring','prm'); if isempty(name), return, end
[list,prms]=gut('getlist','list');
list={list{:},name};
prms={prms{:},gut('getudata','prm')};
[list,i]=unique(list); prms=prms(i);
gut('setlist','list',list,prms)
set(gut('geth','list'),'value',find(strcmp(list,name)));

% callback on delete -----------------------------------------------------------
function cbDel
[list,prms]=gut('getlist','list');
i=gut('getval','list');
list(i)=[]; prms(i)=[];
gut('setlist','list',list,prms)

% callback on clik -------------------------------------------------------------
function cbPrmClk
if strcmp(get(gcf,'selectiontype'),'open') % double click
    set(gcf,'userdata',1);
else
    list=(gut('getstring','list'));
    val=min(gut('getval','list'));
    if ~isempty(list)&~isempty(val), gut('setstring','prm',list{val}); end
end
