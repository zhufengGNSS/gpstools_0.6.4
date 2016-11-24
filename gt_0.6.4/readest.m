function [t,xs,vs,prm,arc]=readest(td,time,type,name,dirs,fb,tunit)
%-------------------------------------------------------------------------------
% [system] : GpsTools
% [module] : read estimation results
% [func]   : read estimation results
% [argin]  : td    = date (mjd-gpst)
%            time  = time vector (sec)
%            type  = estimation type
%            name  = satellite/station name
%           (dirs) = estimation data directory (default:current)
%           (fb)   = estimation direction
%                    ('f':forward,'b':backward,'fb':smoothed)
%           (tunit)= processing unit time (hr) (default:24)
% [argout] : t     = time vector (sec)
%            xs    = estimation results
%            vs    = estimation variences
%            prm   = estimation parameters
%            arc   = arc information
% [note]   :
% [version]: $Revision: 12 $ $Date: 2008-11-25 10:02:15 +0900 (火, 25 11 2008) $
%            Copyright(c) 2004-2006 by T.Takasu, all rights reserved
% [history]: 04/11/17  0.1  new
%            04/12/08  0.2  argin tunit in hr instead of sec
%            05/07/10  0.3  add smoothed estimation
%            06/02/10  0.4  add argout arc
%            06/06/24  0.5  output warning with gt_log()
%-------------------------------------------------------------------------------
if nargin<5, dirs=''; end
if nargin<6, fb='f'; end
if nargin<7, tunit=24; end
t=[]; xs=[]; vs=[]; prm=[]; arc=[]; if isempty(time), return, end
switch fb
case {'f','b'}, [t,xs,vs,vs2,prm,arc]=ReadEstData(td,time,type,name,dirs,fb,tunit);
case 'fb'
    [tf,xsf,vsf,vsf2,prmf,arc]=ReadEstData(td,time,type,name,dirs,'f',tunit);
    [tb,xsb,vsb,vsb2,prmb,arc]=ReadEstData(td,time,type,name,dirs,'b',tunit);
    if length(tf)~=length(tb)|any(tf~=tb)
        gt_log('forward/backward ests incompatible : %s%s',type,fb);
        return
    end
    t=tf; prm=prmf;
    vsf(vsf==0)=nan; vsb(vsb==0)=nan;
    if isempty(vsf2)|isempty(vsb2)
        vs=1./(1./vsf+1./vsb);
        xs=(xsf./vsf+xsb./vsb).*vs;
    else
        xs=repmat(nan,size(xsf)); vs=xs;
        i=1:3; j=4:size(vsf,2);
        for n=1:size(xs,1)
            pf=diag(vsf(n,i))+[0,vsf2(n,1:2);vsf2(n,1),0,vsf2(n,3);vsf2(n,2:3),0];
            pb=diag(vsb(n,i))+[0,vsb2(n,1:2);vsb2(n,1),0,vsb2(n,3);vsb2(n,2:3),0];
            pf=inv(pf); pb=inv(pb); p=inv(pf+pb);
            xs(n,i)=(p*(pf*xsf(n,i)'+pb*xsb(n,i)'))';
            vs(n,i)=diag(p)';
        end
        vs(:,j)=1./(1./vsf(:,j)+1./vsb(:,j));
        xs(:,j)=(xsf(:,j)./vsf(:,j)+xsb(:,j)./vsb(:,j)).*vs(:,j);
    end
otherwise
    gt_log('estimation source error : %s%s',type,fb);
end

% read estimation results data -------------------------------------------------
function [t,xs,vs,vs2,prm,arcs]=ReadEstData(td,time,type,name,dirs,fb,tunit)
t=[]; xs=[]; vs=[]; vs2=[]; prm=[]; arcs=[]; tu=tunit*3600;
for ts=floor(time(1)/tu)*tu:tu:time(end)
    ep=mjdtocal(td,ts);
    file=sprintf('%s%s_%s_%04d%02d%02d%02d.mat',type,fb,name,ep(1:4));
    file=gfilepath(dirs,file,ep,name);
    if exist(file)
        epoch=[]; time=[]; data=[]; covs=[]; covs2=[]; arc=[];
        load(file)
        if ~isempty(epoch)
            [tdd,tss]=caltomjd(epoch);
            time=time+(tdd-td)*86400+tss;
            i=find(ts<=time&time<ts+tu);
            t=[t;time(i)];
            xs=[xs;data(i,:)];
            if ~isempty(covs)
                vs=[vs;covs(i,:)];
            else
                vs=[vs;zeros(length(i),size(data,2))];
            end
            if ~isempty(covs2)
                vs2=[vs2;covs2(i,:)];
            end
            if ~isempty(arc)
                arc(:,1:2)=arc(:,1:2)+(tdd-td)*86400+tss;
                arcs=[arcs;arc];
            end
        end
    else
        gt_log('no estimation results   : %s',file);
    end
end
