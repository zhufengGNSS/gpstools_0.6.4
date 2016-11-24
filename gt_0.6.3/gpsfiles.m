function h=gpsfiles(varargin)
%-------------------------------------------------------------------------------
% [system] : GpsTools
% [module] : filer
% [func]   : filer for gps products/data/estimation results
% [argin]  : none
% [argout] : h   = figure handle
% [note]   :
% [version]: $Revision: 28 $ $Date: 06/07/21 19:09 $
%            Copyright(c) 2004-2006 by T.Takasu, all rights reserved
% [history]: 04/11/15  0.1  new
%            05/06/28  0.2  add funciton to sort files by extentions
%            06/03/06  0.3  add filter, directory up/down buttons
%-------------------------------------------------------------------------------
if nargin<1, h=GpsFilesMain; else feval(varargin{:}); end

% filer main -------------------------------------------------------------------
function h=GpsFilesMain
types={'All(*.*)','Obs Data(obs*.mat)','Sta Position(pos*.mat)','Sat Orbit(eph*.mat)',...
       'Sat/Sta Clock(clk*.mat)','Tropos(zpd*.mat)','ERP(erp*.mat)',...
       'Geocenter(eco*.mat)','Phase Bias(bcp*.mat)','Residuals(res*.mat)',...
       'Multipath Profile(mpp*.mat)','RINEX Obs(*.*o)','RINEX Nav(*.*n)',...
       'SINEX Pos(*.snx)','SP3 Ephem(*.sp3)','SP3 Ephem(*.eph)','RINEX Clk(*.clk)',...
       'SINEX Trop(*.zpd)','IONEX Ion(*.*i)','Figure(*.fig)','Parameter(prm*.m*)',...
       'Log(*.log)','Text(*.txt)'};
filts={'*%p*%t*.*','obs*%p*%t.mat','pos*%p*%t.mat','eph*%p*%t.mat',...
       'clk*%p*%t.mat','zpd*%p*%t.mat','erp*%t.mat',...
       'eco*%t.mat','bcp*%p*%t.mat','res*%t.mat','mpp*%p*%t*.mat'...
       '*%p%t*.*o','*%p%t*.*n','*%p%t*.snx','*%p%t*.sp3','*%p%t*.eph',...
       '*%p%t*.clk','*%p%t*.zpd','*%p%t*.*i','*%p*.fig','prm*.m*','*%p*.log','*%p*.txt'};
prm=loadprm('prm_gpsfiles','prm_gpsfiles_def');
h=gut('newfig','gtfile','GpsTools',[598,400,0,-50],[mfilename,' cbClose'],...
      'resize','on','resizefcn',[mfilename,' cbResize']);
if isempty(h), return, end
hd=gut('newlist','dirs',[0,0,100,370],{},{},1,[mfilename,' cbDirClk' ]);
hf=gut('newlistf','files',[0,0,300,390],{},{},2,[mfilename,' cbFileClk']);
h0=gut('newbtn','rbtn',[0,0,16,16],'\', [mfilename,' cbDirRoot']);
ha=gut('newbtn','abtn',[0,0,16,16],'<', [mfilename,' cbDirUp']);
hb=gut('newbtn','bbtn',[0,0,16,16],'>', [mfilename,' cbDirDown']);
h1=gut('newbtn','dbtn',[0,0,100,16],'Directory', [mfilename,' cbDirGo']);
hc=gut('newbtn','cbtn',[0,0,100,16],'F', [mfilename,' cbFilter']);
h2=gut('newbtn','sbtn1',[0,0,100,16],'Name', [mfilename,' cbFileSort']);
h3=gut('newbtn','sbtn2',[0,0,100,16],'Ext', [mfilename,' cbFileSort']);
h4=gut('newbtn','sbtn3',[0,0,100,16],'Size', [mfilename,' cbFileSort']);
h5=gut('newbtn','sbtn4',[0,0,100,16],'Date/Time', [mfilename,' cbFileSort']);
h6=gut('newbtn','sbtn5',[0,0,100,16],'.', [mfilename,' cbUpdate']);
set(h0,'fontsize',8); set(h1,'fontsize',8); set(h2,'fontsize',8);
set(h3,'fontsize',8); set(h4,'fontsize',8); set(h5,'fontsize',8);
set(h6,'fontsize',8); set(ha,'fontsize',8); set(hb,'fontsize',8);
set(hc,'fontsize',8);
set(gut('newfrm','sep',[0,0,1,1]),'visible','off');
gut('newmenu','file','&File ',...
    {'&New Directory...','-','&Delete ','-','&Uncompact','Com&pact','-','&Close '},{},...
    {[mfilename,' cbNewDir'],'',[mfilename,' cbFileDelete'],'',...
     [mfilename,' cbFileUnzip'],[mfilename,' cbFileZip'],'',[mfilename,' cbClose']});
gut('newmenu','edit','&Edit ',...
    {'&Copy...','&Move...','&Rename...','-','Select &All '},{},...
    {[mfilename,' cbFileCopy'],[mfilename,' cbFileMove'],[mfilename,' cbFileRename'],...
     '',[mfilename,' cbSelAll']});
gut('newmenu','filt','Fi&lter',...
    {'&Filter...','-',types{:}},{'','',filts{:}},[mfilename,' cbFilter']);
gut('newmenu','view','&View ',...
    {'&Update','-','&View...'},{},...
    {[mfilename,' cbUpdate'],'',[mfilename,' cbFileView']});
gut('newmenu','tool','&Tools ',...
    {'&Download...','-','&Options...'},{},...
    {[mfilename,' cbDownload'],'',[mfilename,' cbOpt']});
gut('newcmenu','drvs',[0,0],gut('getdrvlist'),{},[mfilename,' cbDirDrv']);
gut('newcmenu','popup',[0,0],...
    {'&New Directory...','-','&Copy...','&Move...','&Rename...','&Delete ','-',...
     '&Uncompact','Com&pact','-','&View...'},{},...
    {[mfilename,' cbNewDir'],'',[mfilename,' cbFileCopy'],[mfilename,' cbFileMove'],...
     [mfilename,' cbFileRename'],[mfilename,' cbFileDelete'],'',...
     [mfilename,' cbFileUnzip'],[mfilename,' cbFileZip'],'',[mfilename,' cbFileView']});
set(hf,'buttondownfcn',[mfilename,' cbFilePopup'])
set(h,'windowbuttonupfcn',[mfilename,' cbSepUp'])
set(h,'windowbuttonmotionfcn',[mfilename,' cbSepMove'])
set(h,'windowbuttondownfcn',[mfilename,' cbSepDown'])
gut('chkmenu','filt','*%p*%t*.*');
prm.filter={'','',''};
prm.patt='*.*';
set(h,'userdata',prm);
cbResize
[p,f,e]=fileparts(prm.dirs);
if ~isempty(p), cbUpdate(fullfile(p,'..')); else cbUpdate; end

% callback on close ------------------------------------------------------------
function cbClose
prm=get(gcf,'userdata');
prm.pos=get(gcf,'position');
[path,file]=fileparts(which(mfilename));
save(fullfile(path,'settings','prm_gpsfiles.mat'),'prm');
closereq

% callback on separator down/move/up -------------------------------------------
function cbSepDown
p=get(gcf,'position'); pc=get(gcf,'currentpoint'); prm=get(gcf,'userdata'); 
ps=p(3)*prm.possep; if pc(1)<ps-3|ps+3<pc(1), return, end
set(gcf,'pointer','right');
set(gut('geth','sep'),'position',[ps,0,2,p(4)],'visible','on');

function cbSepMove
if ~strcmp(get(gcf,'pointer'),'right'), return, end
p=get(gcf,'position'); pc=get(gcf,'currentpoint');
set(gut('geth','sep'),'position',[pc(1),0,2,p(4)]);

function cbSepUp
if ~strcmp(get(gcf,'pointer'),'right'), return, end
p=get(gcf,'position'); pc=get(gcf,'currentpoint'); prm=get(gcf,'userdata'); 
prm.possep=min(max(pc(1),100)/p(3),0.5);
set(gcf,'pointer','arrow','userdata',prm);
gut('setinv','sep');
cbResize; cbUpdate;

% callback on resize -----------------------------------------------------------
function cbResize
prm=get(gcf,'userdata'); if isempty(prm), return, end
p=get(gcf,'position');
w=max(p(3)+2,20); hh=max(p(4)-16,10);
set(gut('geth','dirs'), 'position',[0,0,w*prm.possep-1,hh]);
set(gut('geth','files'),'position',[w*prm.possep+1,0,max(w*(1-prm.possep)-1,10),hh]);
set(gut('geth','rbtn'),'position',[2,hh,16,16]);
set(gut('geth','drvs'),'position',[2,hh]);
set(gut('geth','abtn'),'position',[18,hh,16,16]);
set(gut('geth','bbtn'),'position',[34,hh,16,16]);
set(gut('geth','dbtn'),'position',[50,hh,max(w*prm.possep-54,10),16]);
w2=60; w3=106; w4=16; w1=max(w*(1-prm.possep)-w2-w3-23,60); 
set(gut('geth','cbtn'),'position',[w*prm.possep+3,hh,16,16]);
set(gut('geth','sbtn1'),'position',[w*prm.possep+19,hh,w1-51,16]);
set(gut('geth','sbtn2'),'position',[w*prm.possep+w1-32,hh,36,16]);
set(gut('geth','sbtn3'),'position',[w*prm.possep+w1+4,hh,w2,16]);
set(gut('geth','sbtn4'),'position',[w*prm.possep+w1+w2+4,hh,w3,16]);
set(gut('geth','sbtn5'),'position',[w*prm.possep+w1+w2+w3+4,hh,w4,16]);

% update directory -------------------------------------------------------------
function cbUpdate(dirs)
h=findobj(0,'tag','gtfile'); if isempty(h), return, end
f=gcf; figure(h);
if nargin<1, dirs=gut('getdirlist','dirs'); end
if (~isempty(dirs)&exist(dirs)~=7&~strcmp(dirs,'\'))|...
   strcmp(dirs,'..')|strcmp(dirs,'.'), return, end
gut('updatedirlist','dirs',dirs);
prm=get(gcf,'userdata'); prm.dirs=dirs; set(gcf,'userdata',prm);
UpdateFiles(0);
figure(f);

% update file ------------------------------------------------------------------
function UpdateFiles(opt)
prm=get(gcf,'userdata');
patt=FilePatt(prm.patt,prm.filter);
while 1, p=strrep(patt,'**','*'); if strcmp(p,patt), break, else patt=p; end, end
if strcmp(patt,'*.*'), s=''; else s=[' (',patt,')']; end
set(gcf,'name',['GpsTools : ',prm.dirs,s]);
[p,d,e]=fileparts(prm.dirs);
if strcmp(d,'.')&strcmp(e,'.')
    gut('setstring','files',{});
else
    siz0=gut('getsize','sbtn1'); siz1=gut('getsize','sbtn2');
    siz2=gut('getsize','sbtn3'); siz3=gut('getsize','sbtn4');
    gut('updatefilelist','files',prm.dirs,patt,prm.sort,prm.ud,...
        [16+siz0(1)+siz1(1)+4,siz2(1)+2,siz3(1)+2],opt)
end
fa={'normal','normal','normal','normal'}; fa{prm.sort+1}='italic';
set(gut('geth','sbtn1'),'fontangle',fa{1});
set(gut('geth','sbtn2'),'fontangle',fa{4});
set(gut('geth','sbtn3'),'fontangle',fa{2});
set(gut('geth','sbtn4'),'fontangle',fa{3});

% callback on button directory root --------------------------------------------
function cbDirRoot
gut('setvis','drvs');

% callback on menu select drive ------------------------------------------------
function cbDirDrv
cbUpdate([get(gcbo,'label'),'..']);

% callback on button directory up ----------------------------------------------
function cbDirUp
dirs=gut('getdirlist','dirs');
[parent,file,ext]=fileparts(dirs);
if ~isempty(parent), cbUpdate(parent); end

% callback on button directory down --------------------------------------------
function cbDirDown
dirs=gut('getdirlist','dirs');
if isempty(dirs)|dirs(end)=='.', return, end
if dirs(end)=='\'; dirs=[dirs,'..']; else dirs=[dirs,'\..']; end
cbUpdate(dirs)

% callback on button directory click -------------------------------------------
function cbDirGo
gut('newdlg','','Go To Directory',[300,54]);
set(gut('newedit','go',[7,30,288,20],'',1),'callback',[mfilename,' cbOk']);
gut('newokcancelbtn','',[115,3,180,23]);
if ~gut('waitok'), return, end
dirs=gut('getstring','go');
close
cbUpdate(dirs)

function cbOk, set(gcf,'userdata',1);

% callback on directory click --------------------------------------------------
function cbDirClk
dirs=gut('getdirlist','dirs');
if strcmp(get(gcf,'selectiontype'),'open')
    [parent,file,ext]=fileparts(dirs);
    if strcmp(file,'.'), dirs=parent;
    elseif dirs(end)=='\'; dirs=[dirs,'..']; else dirs=[dirs,'\..']; end
end
cbUpdate(dirs)

% callback on file popupmenu ----------------------------------------------------
function cbFilePopup
if ~strcmp(get(gcf,'selectiontype'),'alt'), return, end
set(gut('geth','popup'),'position',get(gcf,'currentpoint'))
gut('showmenu','popup');

% callback on file click -------------------------------------------------------
function cbFileClk
if ~strcmp(get(gcf,'selectiontype'),'open'), return, end
[dirs,isdir]=gut('getfilelist','files','');
if ~isempty(dirs)&isdir(1)
    parent=gut('getdirlist','dirs');
    if ~strcmp(dirs{1},'..')
        cbUpdate(fullfile(parent,dirs{1}));
        set(gcbo,'value',find(strcmp(get(gcbo,'string'),'<..>')));
    else
        [parent,f,e]=fileparts(parent);
        cbUpdate(parent);
        set(gcbo,'value',find(strcmp(get(gcbo,'string'),['<',f,e,'>'])));
    end
else, cbFileView; end

% callback on file sort --------------------------------------------------------
function cbFileSort
prm=get(gcf,'userdata');
switch get(gcbo,'tag')
case 'sbtn1', sort=0;
case 'sbtn2', sort=3;
case 'sbtn3', sort=1;
case 'sbtn4', sort=2; end
if sort~=prm.sort, prm.sort=sort; elseif prm.ud, prm.ud=0; else prm.ud=1; end
set(gcf,'userdata',prm);
UpdateFiles(1);

% callback on filter -----------------------------------------------------------
function cbFilter
prm=get(gcf,'userdata');
patt=get(gcbo,'userdata');
if isempty(patt)
    prm1={
    'e','Pattern',      prm.filter{1}
    'e','Sats/Stations',prm.filter{2}
    'y','Date',         prm.filter{3}
    };
    gut('newdlg','','Filter',[215,105])
    gut('newprms','prm1',[10,32,200,23,116],prm1);
    gut('newokcancelbtn','',[50,4,160,22]);
    if ~gut('waitok'), return, end
    [prm.filter{1},prm.filter{2},prm.filter{3}]=gut('getprms','prm1');
    close
else
    gut('unchkmenu','filt');
    set(gcbo,'checked','on');
    prm.patt=patt;
end
set(gcf,'userdata',prm);
cbUpdate;

% file pattern -----------------------------------------------------------------
function patt=FilePatt(patt,filter)
if ~isempty(filter{1}), patt=['*',filter{1},'*']; end
if length(filter{3})==3
    td=caltomjd(filter{3});
    switch patt
    case {'*%p%t*.*o','*%p%t*.*n','*%p%t*.*i'}
        date=sprintf('*%03d*',td-caltomjd([filter{3}(1),1,1])+1);
    case {'*%p%t*.sp3','*%p%t*.eph','*%p%t*.clk','*%p%t*.zpd'}
        d=td-44244; w=floor(d/7); date=sprintf('*%04d%d*',w,floor(d-w*7));
    case '*%p%t*.snx'
        d=td-44244; w=floor(d/7); date=sprintf('*%02dP%04d*',mod(filter{3}(1),100),w);
    otherwise date=sprintf('*%04d%02d%02d*',filter{3});
    end
else date=''; end
patt=strrep(patt,'%t',date);
patt=strrep(patt,'%p',filter{2});

% callback on select all -------------------------------------------------------
function cbSelAll
h=gut('geth','files');
files=get(h,'string');
sel=ones(length(files),1);
sel(strmatch('<',files))=0;
set(h,'value',find(sel));

% callback on new directory ----------------------------------------------------
function cbNewDir
ndir=gut('newdirdlg'); if isempty(ndir), return, end
dirs=gut('getdirlist','dirs');
wd=pwd; cd(dirs),
if mkdir(ndir)==1, cbUpdate, end
cd(wd)

% callback on copy -------------------------------------------------------------
function cbFileCopy
files=gut('getfilelist','files','');
if isempty(files), return, end
sdir=gut('getdirlist','dirs');
[ok,ddir]=gut('dirdlg',sdir); if ~ok, return, end
gut('newmsgbar','gfmsg','',[480,50],[0,0.7,1],1);
for n=1:length(files)
    if files{n}(1)~='<'
        [stat,msg]=dos(['copy "',fullfile(sdir,files{n}),'" "',ddir,'"']);
        if stat~=0&~gut('confdlg',['Error : ',msg],{'OK','Abort'}), break, end
        if gut('setmsgbar','gfmsg',['Copying : ',files{n}],n/length(files)), break, end
    end
end
gut('closemsgbar','gfmsg');
UpdateFiles(1);

% callback on options ----------------------------------------------------------
function cbOpt
sel1={{'.Z','.gz'},0:1};
prm=get(gcf,'userdata');
prm1={
'p','Compression Suffix',sel1
};
prm2={
' ','Text Viewer', '';'f','',''
' ','Image Viewer','';'f','',''
' ','Browser',     '';'f','',''
};
gut('newdlg','','Options',[350,178]);
gut('newprms','prm1',[8,150,337,24,90],prm1);
gut('newprms','prm2',[8,34,337,18,90],prm2);
gut('newokcancelbtn','',[165,4,180,23]);
gut('setprms','prm1',prm.compsufx);
gut('setprms','prm2','',prm.cmd1,'',prm.cmd2,'',prm.cmd3);
if ~gut('waitok'), return, end
[prm.compsufx]=gut('getprms','prm1');
[q,prm.cmd1,q,prm.cmd2,q,prm.cmd3]=gut('getprms','prm2');
close
set(gcf,'userdata',prm);
cbResize

% callback on move -------------------------------------------------------------
function cbFileMove
[files,isdir]=gut('getfilelist','files','');
if isempty(files), return, end
sdir=gut('getdirlist','dirs');
[ok,ddir]=gut('dirdlg',sdir); if ~ok, return, end
gut('newmsgbar','gfmsg','',[480,50],[0,0.7,1],1);
for n=1:length(files)
    if ~isdir(n)
        [stat,msg]=dos(['move "',fullfile(sdir,files{n}),'" "',ddir,'"']);
        if stat~=0&~gut('confdlg',['Error : ',files{n},' ',msg],{'OK','Abort'}), break, end
        if gut('setmsgbar','gfmsg',['Moving : ',files{n}],n/length(files)), break, end
    end
end
gut('closemsgbar','gfmsg');
cbUpdate

% callback on rename -----------------------------------------------------------
function cbFileRename
[files,isdir]=gut('getfilelist','files','');
if isempty(files)|isdir(1), return, end
sdir=gut('getdirlist','dirs');
gut('newdlg','','Rename',[240,55]);
h=gut('newedit','name',[5,30,230,20],files{1},1);
set(h,'callback',[mfilename,' cbFileRenameOk']);
gut('newokcancelbtn','',[75,5,160,20]);
if ~gut('waitok'), return, end
newfile=gut('getstring','name');
close
[stat,msg]=dos(['move "',fullfile(sdir,files{1}),'" "',fullfile(sdir,newfile),'"']);
if stat==0, cbUpdate; else gut('errdlg',['Error : ',files{1},' ',msg]); end

function cbFileRenameOk, set(gcf,'userdata',1);

% callback on delete -----------------------------------------------------------
function cbFileDelete
[files,isdir]=gut('getfilelist','files','');
if isempty(files)|~gut('confdlg','Really Delete Files/Directories?'), return, end
sdir=gut('getdirlist','dirs');
gut('newmsgbar','gfmsg','',[480,50],[0,0.7,1],1);
for n=1:length(files)
    if ~isdir(n), delete(fullfile(sdir,files{n}));
    elseif ~strcmp(files{n},'..')
        [ok,msg]=rmdir(fullfile(sdir,files{n}));
        if ~ok&~gut('confdlg',['Error : ',files{n},' ',msg],{'OK','Abort'}), break, end
    end
    if gut('setmsgbar','gfmsg',['Deleting : ',files{n}],n/length(files)), break, end
end
gut('closemsgbar','gfmsg');
cbUpdate

% callback on unzip/uncompact --------------------------------------------------
function cbFileUnzip
[files,isdir]=gut('getfilelist','files','');
if isempty(files), return, end
dirs=gut('getdirlist','dirs');
gut('newmsgbar','gfmsg','',[480,50],[0,0.7,1],1);
for n=1:length(files)
    if ~isdir(n)
        [f,stat,msg]=uncompact(fullfile(dirs,files{n}));
        if stat~=0&~gut('confdlg',['Error : ',files{n},' ',msg],{'OK','Abort'}), break, end
        if gut('setmsgbar','gfmsg',['Uncompacting : ',files{n}],n/length(files)), break, end
    end
end
gut('closemsgbar','gfmsg');
cbUpdate

% callback on zip/compact ------------------------------------------------------
function cbFileZip
prm=get(gcf,'userdata'); if prm.compsufx, opt='gz'; else opt=''; end
[files,isdir]=gut('getfilelist','files','');
if isempty(files), return, end
dirs=gut('getdirlist','dirs');
gut('newmsgbar','gfmsg','',[480,50],[0,0.7,1],1);
for n=1:length(files)
    if ~isdir(n)
        [f,stat,msg]=compact(fullfile(dirs,files{n}),opt);
        if stat~=0&~gut('confdlg',['Error : ',files{n},' ',msg],{'OK','Abort'}), break, end
        if gut('setmsgbar','gfmsg',['Compacting : ',files{n}],n/length(files)), break, end
    end
end
gut('closemsgbar','gfmsg');
cbUpdate

% callback on download ---------------------------------------------------------
function cbDownload
dirs=gut('getdirlist','dirs');
gpsdown; gpsdown('cbUpdate',dirs);

% callback on view -------------------------------------------------------------
function cbFileView
files=gut('getfilelist','files',gut('getdirlist','dirs'));
if isempty(files), return, end
prm=get(gcf,'userdata');
[path,file,ext]=fileparts(files{1});
switch lower(ext)
case '.mat',
    if any(strncmp(file,{'mpc_','mpp_'},2))
        plotmp(files{1});
        return;
    end
    str_={' Name                 Size          Bytes  Class   Value',''};
    file_=files{1};
    for s_=whos('-file',file_)'
        load(files{1},s_.name);
        ss_=ViewVar(s_.name,eval(s_.name),'');
        str_={str_{:},ss_{:}};
    end
    ViewerS(file_,str_)
case '.fig'
    openfig(files{1});
case {'.jpg','.jpeg','.gif','.tif','.tiff','.bmp','.png','.emf','.wmf'}
    if ~isempty(prm.cmd2)
        [s,w]=dos(['"',prm.cmd2,'" "',files{1},'" &']);
        if s, gut('errdlg',w); end
    else
        try, [img,map]=imread(files{1}); catch gut('errdlg',lasterr); img=[]; end
        if ~isempty(img)
            w=size(img,2); h=size(img,1);
            gut('newfig','',files{1},[w,h]);
            if ~isempty(map), set(gcf,'colormap',map); end
            image(img); axis off;
            p=get(gcf,'position');
            set(gca,'position',[(1-w/p(3))/2,(1-h/p(4))/2,w/p(3),h/p(4)]);
        end
    end
case {'.html','.htm'}
    if ~isempty(prm.cmd3)
        [s,w]=dos(['"',prm.cmd3,'" "',files{1},'" &']);
        if s, gut('errdlg',w); end
    end
case {'.exe','.dll','.zip'}
    ;
otherwise
    if isempty(prm.cmd1), Viewer(files);
    else dos(['"',prm.cmd1,'" "',files{1},'" &']); end
end

% view variable ----------------------------------------------------------------
function ss=ViewVar(name,value,head)
s=whos('value');
siz=sprintf('%d',s.size(1));
for n=2:length(s.size), siz=sprintf('%sx%d',siz,s.size(n)); end
ss={sprintf(' %-20s %-8s %10d  %-7s',[head,name],siz,s.bytes,s.class)};
switch s.class
case 'double', ss={[ss{:},' ','[',NumToStr(value),']']};
case 'cell',   ss={[ss{:},' ','{',CellToStr(value),'}']};
case 'char',   ss={[ss{:},' ''',value,'''']};
case 'struct'
    for f=fieldnames(value)'
        sss=ViewVar(['.',f{1}],eval(['value(1).',f{1}]),[head,' ']);
        ss={ss{:},sss{:}};
    end
    if length(value)>1, ss={ss{:},[head,'  ...']}; end
end

function str=ToStr(value)
if ischar(value), str=['''',value,''''];
elseif isnumeric(value), str=['[',NumToStr(value),']'];
else str='...'; end

function str=NumToStr(value)
siz=size(value); str='';
for n=1:siz(1)
    if n>1, str=[str,';']; end
    for m=1:siz(2)
        if m>1, str=[str,',']; end
        if length(str)>120, str=[str,'...']; return, end
        str=[str,num2str(value(n,m))];
    end
end

function str=CellToStr(value)
siz=size(value); str='';
for n=1:siz(1)
    if n>1, str=[str,';']; end
    for m=1:siz(2)
        if m>1, str=[str,',']; end
        if length(str)>120, str=[str,'...']; return, end
        str=[str,ToStr(value{n,m})];
    end
end

% show by text viewer ----------------------------------------------------------
function Viewer(files)
figure(gut('newviewer','','Text Viewer ',[600,500],'file',files));

function ViewerS(title,str)
figure(gut('newviewer','',title,[600,500],'str',str));
