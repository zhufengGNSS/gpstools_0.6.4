function t_readeph

td=caltomjd([2003,12,1]);
time=0:900:86400-900; % GPST
sat='GPS01';

erp=ReadErp(td,time-13,-32,'erp');
U=zeros(3,3,length(time)); dUdt=U;
for n=1:length(time)
    [U(:,:,n),a,b,c,d,e,dUdt(:,:,n)]=EcsfToEcef(td+(time(n)-13)/86400,erp(n,:),-32);
end
ephs=ReadEph(td,time,sat,U,'eph','nima'); % ‰q¯ˆÊ’u‘¬“x

[ep,ti,ephe]=ReadSp3('eph\nim12471.eph');
ephr=zeros(length(time),6);
for n=1:length(time) % ECEF->ECI
    ephr(n,:)=[ephe(n,1:3,1)*U(:,:,n),...
               ephe(n,5:7,1)*U(:,:,n)+ephe(n,1:3,1)*dUdt(:,:,n)];
end

subplot(2,1,1), plot(time/3600,ephs(:,1:3)-ephr(:,1:3))
subplot(2,1,2), plot(time/3600,ephs(:,4:6)-ephr(:,4:6))
