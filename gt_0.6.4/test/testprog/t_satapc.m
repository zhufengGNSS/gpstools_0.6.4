function t_satapc
% test for satapc

[apc1,apc2,apv1,apv2,stat]=readpcv('','GPS23','..\..\data\igs05.atx');

apc1=apc1';
apc2=apc2';
%ecc=[0;0;0];
%apc1=[0;0;0];
%apc2=[0;0;0];
apv1=shiftdim(apv1);
apv2=shiftdim(apv2);

az=(0:36:360)*pi/180;
nadir=(0:1:16)*pi/180;

for i=1:length(az)
    msg=sprintf('%3.0f: ',az(i)*180/pi);
    for j=1:length(el)
        apr=satapc([az(i),nadir(j)],apc1,apc2,apv1,apv2);
        msg=[msg,sprintf(' %7.4f',apr(1))];
    end
    disp(msg);
end
for i=1:length(az)
    msg=sprintf('%3.0f: ',az(i)*180/pi);
    for j=1:length(el)
        apr=satapc([az(i),nadir(j)],apc1,apc2,apv1,apv2);
        msg=[msg,sprintf(' %7.4f',apr(2))];
    end
    disp(msg);
end
