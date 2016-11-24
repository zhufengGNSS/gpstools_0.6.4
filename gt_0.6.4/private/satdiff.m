function [dzs,G,R,iz]=satdiff(dz,dgds,sig,ig,ircvs)
%-------------------------------------------------------------------------------
% [system] : GpsTools
% [module] : single difference between satellites
% [func]   : generate single difference measurement between satellites
% [argin]  : dz  = residuals
%            dgds= partial derivatives of model
%            sig = meas. noise std. deviations
%            ig  = meas. index [t,isat,ircv,arcf;...]
%            ircvs = baselines of stations [ircv1,ircv2,...] (no use)
% [argout] : dzs = difference of residuals
%            G   = difference of partial derivatives
%            R   = meas. noise covariences
%            iz  = meas. index [isat1,isat2,ircv,0;...]
% [note]   :
% [version]: $Revision: 12 $ $Date: 2008-11-25 10:02:15 +0900 (火, 25 11 2008) $
% [history]: 04/11/19  0.1  new
%-------------------------------------------------------------------------------
nz=0; dzs=zeros(1000,1); G=zeros(1000,size(dgds,2)); R=zeros(1000,1000);
iz=zeros(1000,4);
nrcv=max(max(ig(:,3)));

for n=1:nrcv
    i=find(ig(:,3)==n);
    for j=1:length(i)-1
        nz=nz+1;
        dzs(nz)=dz(i(j))-dz(i(j+1));
        G(nz,:)=dgds(i(j),:)-dgds(i(j+1),:);
        R(nz,nz)=sig(i(j))^2+sig(i(j+1))^2;
        if j>1, R(nz-1,nz)=sig(i(j))^2; R(nz,nz-1)=sig(i(j))^2; end
        iz(nz,:)=[ig(i(j),2),ig(i(j+1),2),n,0];
    end
end
dzs=dzs(1:nz); G=G(1:nz,:); R=R(1:nz,1:nz); iz=iz(1:nz,:);
