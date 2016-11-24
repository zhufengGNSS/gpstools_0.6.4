function makemex(varargin)
%-------------------------------------------------------------------------------
% [system] : GpsTools
% [module] : mex make
% [func]   : make mex routines
% [argin]  :(targets) = target mex routines
% [argout] : none
% [note]   :
% [version]: $Revision: 14 $ $Date: 06/07/08 14:12 $
%            Copyright(c) 2004-2006 by T.Takasu, all rights reserved
% [history]: 04/11/19  0.1  new
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
{'filterekf.c'}
{'filtersrcf.c'}
{'geodtoecef.c','geodtoecef_e.c'}
{'ion_klob.c'}
{'isrange.c','isrange_e.c'}
{'lambda.c'}
{'mapf_nmf.c','mapf_nmf_e.c','caltomjd_e.c','mjdtocal_e.c'}
{'mapf_gmf.c','gmf.obj'}
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
 'ecsftosatf_e.c','gmf.obj'}
{'readgrib.c'}
{'readdgrb.c'}
{'readrinexclk.c','caltomjd_e.c'}
{'readrinexnav.c'}
{'readrinexobs.c','caltomjd_e.c'}
{'readsp3.c','caltomjd_e.c'}
{'readgsipos.c','caltomjd_e.c'}
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
};
if nargin==0, index=1:length(srcs);
else
    index=[]; for n=1:length(srcs)
        if any(strcmp(srcs{n}{1}(1:end-2),varargin)), index=[index,n]; end
    end
end
switch lower(computer)
case 'pcwin',  ext='.dll';    opts='';     % windows
case 'glnx86', ext='.mexglx'; opts='';     % linux
end
for n=index
    obj=dir([srcs{n}{1}(1:end-2),ext]); d='';
    cmd=['mex ',opts,' ']; make=0;
    for m=1:length(srcs{n})
        src=dir(['mexsrc',filesep,srcs{n}{m}]);
        if ~isempty(src)
            cmd=[cmd,' mexsrc',filesep,srcs{n}{m}];
            if isempty(obj)|datenum(obj.date)<datenum(src.date), make=1; end
        else
            disp([srcs{n}{m},' no src']), break;
        end
    end
    if make|isempty(obj)
        disp(cmd), eval(cmd)
    else
        disp([obj.name,' up to date']);
    end
end
