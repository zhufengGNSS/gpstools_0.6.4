function t_rcvapc
% test for rcvapc

[apc1,apc2,apv1,apv2,stat]=readpcv('','ASH701933B_M','..\..\data\igs05.atx');

ecc=[0.9;0.6;0.5];
apc1=apc1';
apc2=apc2';
%ecc=[0;0;0];
%apc1=[0;0;0];
%apc2=[0;0;0];
apv1=shiftdim(apv1);
apv2=shiftdim(apv2);

az=(0:36:360)*pi/180;
el=(0:6:90)*pi/180;

for i=1:length(az)
    msg=sprintf('%3.0f: ',az(i)*180/pi);
    for j=1:length(el)
        apr=rcvapc([az(i),el(j)],apc1,apc2,ecc,apv1,apv2);
        msg=[msg,sprintf(' %7.4f',apr(1))];
    end
    disp(msg);
end
for i=1:length(az)
    msg=sprintf('%3.0f: ',az(i)*180/pi);
    for j=1:length(el)
        apr=rcvapc([az(i),el(j)],apc1,apc2,ecc,apv1,apv2);
        msg=[msg,sprintf(' %7.4f',apr(2))];
    end
    disp(msg);
end
