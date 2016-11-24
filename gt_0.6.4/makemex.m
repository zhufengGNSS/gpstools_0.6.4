function makemex(varargin)
%-------------------------------------------------------------------------------
% [system] : GpsTools
% [module] : mex make
% [func]   : make mex routines
% [argin]  : (options) = -outdir dir: specify output directory
%                        -v         : display detaied log
%                        -all       : re-generate all
%            (targets) = target mex routines
% [argout] : none
% [note]   :
% [version]: $Revision: 20 $ $Date: 2009-05-01 04:15:33 +0900 (金, 01 5 2009) $
%            Copyright(c) 2004-2008 by T.Takasu, all rights reserved
% [history]: 04/11/19  0.1  new
%            08/11/21  0.2  support matlab 7, extension .dll->mexw32 (gt_0.6.4)
%                           add trop_gpt readvmfgrid mapf_vmf1
%                           delete filterekf,filtersrcf,lambda
%                           support windows 64bit
%-------------------------------------------------------------------------------
srcs={
{'caltomjd.c','caltomjd_e.c'}
{'dbldiff.c','dbldiff_e.c'}
{'eceftogeod.c','eceftogeod_e.c'}
{'ecsftoecef.c','ecsftoecef_e.c','ceppred.obj'}
{'ecsftosatf.c','ecsftosatf_e.c'}
{'eletostate.c','eletostate_e.c'}
{'ephpl.c','ephpl_e.c'}
{'satfixed.c','satfixed_e.c'}
{'geodtoecef.c','geodtoecef_e.c'}
{'ion_klob.c'}
{'isrange.c','isrange_e.c'}
{'mapf_nmf.c','mapf_nmf_e.c','caltomjd_e.c','mjdtocal_e.c'}
{'mapf_gmf.c','gmf.obj'}
{'mapf_vmf1.c','vmf1_ht.obj'}
{'mjdtocal.c','mjdtocal_e.c'}
{'navtostate.c','navtostate_e.c','caltomjd_e.c'}
{'ntrp113.c'}
{'phwindup.c','phwindup_e.c','eceftogeod_e.c','geodtoecef_e.c'}
{'pointpos.c','pointpos_e.c','navtostate_e.c','satazel_e.c','geodtoecef_e.c','eceftogeod_e.c',...
 'trop_saast_e.c','caltomjd_e.c','qtcmn.c'}
{'rcvapc.c','rcvapc_e.c'}
{'rcvmpc.c','rcvmpc_e.c','sphfunc_e.c'}
{'rangemodel.c','rangemodel_e.c','eceftogeod_e.c','satrange_e.c','satazel_e.c',...
 'trop_saast_e.c','mapf_nmf_e.c','satapc_e.c','rcvapc_e.c','relcorr_e.c','phwindup_e.c',...
 'rcvmpc_e.c','geodtoecef_e.c','mjdtocal_e.c','caltomjd_e.c','sphfunc_e.c',...
 'ecsftosatf_e.c','gmf.obj','vmf1_ht.obj'}
{'readgrib.c'}
{'readdgrb.c'}
{'readrinexclk.c','caltomjd_e.c'}
{'readrinexnav.c'}
{'readrinexobs.c','caltomjd_e.c'}
{'readsp3.c','caltomjd_e.c'}
{'readgsipos.c','caltomjd_e.c'}
{'readvmfgrid.c'}
{'relcorr.c','relcorr_e.c'}
{'satorbit.c','satorbit_e.c','ecsftoecef_e.c','ephpl_e.c','geograv_e.c','shadowfunc_e.c','ceppred.obj'}
{'satorbit_s.c','satorbit_e.c','ecsftoecef_e.c','ephpl_e.c','geograv_e.c','shadowfunc_e.c','ceppred.obj'}
{'satrange.c','satrange_e.c'}
{'satapc.c','satapc_e.c'}
{'satazel.c','satazel_e.c','eceftogeod_e.c','geodtoecef_e.c'}
{'geodist.c','navtostate_e.c','caltomjd_e.c'}
{'shadowfunc.c','shadowfunc_e.c'}
{'sitedisp.c','sitedisp_e.c'}
{'sphfunc.c','sphfunc_e.c'}
{'state_2body.c'}
{'state_geoj2.c','ecsftoecef_e.c','ceppred.obj'}
{'state_kepler.c','eletostate_e.c','statetoele_e.c','qtcmn.c'}
{'state_rcvclock.c'}
{'state_satclock.c'}
{'statetoele.c','statetoele_e.c'}
{'trop_saast.c','trop_saast_e.c'}
{'trop_gpt.c','gpt.obj'}
};
n=1; opts=''; all=0;
while n<=nargin
    if any(strcmp(varargin{n},{'-f','-outdir'}))&i<nargin
        opts=[opts,' ',varargin{n},' ',varargin{n+1}]; n=n+1;
    elseif strcmp(varargin{n},'-all')
        all=1;
    elseif strncmp(varargin{n},'-',1)
        opts=[opts,' ',varargin{n}];
    end
    n=n+1;
end

switch lower(computer)
case 'pcwin',  eobj='_w32'; ext='.mexw32'; cfg='mexopts_w32.bat'; % windows 32bit
case 'pcwin64',eobj='_w64'; ext='.mexw64'; cfg='mexopts_w64.bat'; % windows 64bit
case 'glnx86', eobj='';     ext='.mexglx'; cfg='mexopts_linux.sh'; % linux
otherwise, disp(['no supported archtecture: ',computer]); return;
end
for n=1:length(srcs)
    [d,f,e]=fileparts(srcs{n}{1}); obj=[d,f,ext];
    cmd=['mex -f ',cfg,' -output ',obj,' ',opts]; make=0;
    sobj=dir(obj);
    for m=1:length(srcs{n})
        [d,f,e]=fileparts(srcs{n}{m});
        if strcmp(e,'.obj'), srcs{n}{m}=[d,f,eobj,e]; end
        ssrc=dir(['mexsrc',filesep,srcs{n}{m}]);
        if ~isempty(ssrc)
            cmd=[cmd,' mexsrc',filesep,srcs{n}{m}];
            if isempty(sobj)|datenum(sobj.date)<datenum(ssrc.date), make=1; end
        else
            disp([srcs{n}{m},' no src']), break;
        end
    end
    if all|make
        try
            disp(cmd), eval(cmd)
        catch
            disp(['cmd exec error ',cmd])
        end
    else
        disp([sobj.name,' up to date']);
    end
end
