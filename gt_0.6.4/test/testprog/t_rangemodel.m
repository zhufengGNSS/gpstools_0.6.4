function t_rangemodel(ver)
% unit test of rangemodel
if nargin<1, ver='0.6.4'; end
i=0; j=0;
load('t_rangemodel_ws0.mat');
whos

if strcmp(ver,'0.5.6')
    ants=[prm.ants(1:3,:);zeros(16,size(prm.ants,2))];
    antp=antp(:,[1:9,10:28,1397:1415]);
    [y_,H_,azel_,phn_,drds_,j_]=   rangemodel(td,t,x,ix,[tz(i),double(iz(i,:))],...
        state.ephs,state.clks'*C,state.clkr'*C,state.posr,state.dpos,...
        state.U,state.dx,state.dy,state.du,state.rsun,ants,antp',...
        phs,prm.elmin,prm.elmax,state.metr,prm.trop,prm.mapf,0,prm.mpcc,prm.mpcs,...
        prm.f1,prm.f2,prm.corrf,0,prm.rattmodel);
else
    [y_,H_,azel_,phn_,drds_,j_,ds_]=rangemodel(td,t,x,ix,[tz(i),double(iz(i,:))],...
        state.ephs,state.clks'*C,state.clkr'*C,state.zpdr,state.posr,state.dpos,...
        state.bcpr,state.U,state.dx,state.dy,state.du,state.rsun,prm.ants,antp',...
        phs,prm.elmin,prm.elmax,state.metr,prm.trop,prm.mapf,prm.mpcc,prm.mpcs,...
        prm.f1,prm.f2,prm.corrf,double(strcmp(prm.rcvs,prm.clkref)),prm.rattmodel,...
        double(prm.rposmodel==3),state.mfcr);
end


load('t_rangemodel_ws1.mat');

% differnce
if strcmp(ver,'0.5.6')
    disp('ver.0.6.4 result=')
    for k=i'
        disp(sprintf('%.5f ',ds(k,3:end)));
    end
else
    ds_-ds(i,3:end)
end
y_-y,
H_-H,
azel_-azel,
phn_-phn,
drds_-drds,
j_-j
