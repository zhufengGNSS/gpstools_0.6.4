function gengeoidegm96(file)
if nargin<1, file='..\..\geoid\egm96\WW15MGH.GRD'; end
[gprm.lat2,gprm.lat1,gprm.lon1,gprm.lon2,gprm.dy,gprm.dx]=...
    textread(file,'%f %f %f %f %f %f',1);
gprm.nx=(gprm.lon2-gprm.lon1)/gprm.dx+1;;
gprm.ny=(gprm.lat2-gprm.lat1)/gprm.dy+1;
x=dlmread(file,' ',1,0);
data=[];
for n=1:181:size(x,1)
    xx=reshape(x(n:n+180,1:8)',1,181*8);
    data=[data;xx(1:gprm.nx)];
end
data=single(flipud(data));
gprm.type =0;
gprm.rcflg=[0,0,0];
gprm.smode=[0,0];
gprm.lov  =0;
gprm.lati1=0;
gprm.lati2=0;
gprm.latp =0;
gprm.lonp =0;
gprm.x0   =1;
gprm.y0   =1;
layer     =0;
save('geoid_egm96.mat','gprm','layer','data');
