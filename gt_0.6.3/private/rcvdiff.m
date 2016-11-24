function [dzs,G,R,iz]=rcvdiff(dz,dgds,sig,ig,ircvs)
%-------------------------------------------------------------------------------
% [system] : GpsTools
% [module] : single difference between stations
% [func]   : generate sigle difference measurement between stations
% [argin]  : dz  = residuals
%            dgds= partial derivatives of model
%            sig = meas. noise. std. deviations
%            ig  = meas. index [t,isat,ircv,arcf;...]
%            ircvs = baselines of stations [ircv1,ircv2;...]
% [argout] : dzs = difference of residuals
%            G   = difference of partial derivatives
%            R   = meas. noise covariences
%            iz  = meas. index [isat,0,ircv1,ircv2;...]
% [note]   :
% [version]: $Revision: 1 $ $Date: 06/03/25 22:04 $
% [history]: 04/11/02  0.2  new
%-------------------------------------------------------------------------------
nz=0; dzs=zeros(1000,1); G=zeros(1000,size(dgds,2)); id=zeros(1000,2);
iz=zeros(1000,4);

for n=1:size(ircvs,1) % for each baseline
    i1=find(ig(:,3)==ircvs(n,1));
    i2=find(ig(:,3)==ircvs(n,2));
    [ii,j1,j2]=intersect(ig(i1,2),ig(i2,2));
    for k=1:length(ii)
        nz=nz+1; i11=i1(j1(k)); i12=i2(j2(k));
        dzs(nz)=dz(i11)-dz(i12);
        G(nz,:)=dgds(i11,:)-dgds(i12,:);
        id(nz,:)=[i11,i12];
        iz(nz,:)=[ii(k),0,ircvs(n,:)];
    end
end
dzs=dzs(1:nz); G=G(1:nz,:); iz=iz(1:nz,:); R=zeros(nz);
for n=1:nz
    i=find(id(n,1)==id(:,1)|id(n,1)==id(:,2));
    j=find(id(n,2)==id(:,1)|id(n,2)==id(:,2));
    R(n,i)=R(n,i)+sig(id(n,1))^2;
    R(n,j)=R(n,j)+sig(id(n,2))^2;
end
