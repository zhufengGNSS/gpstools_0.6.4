function t_readantex
% test for readantex

[ants,apc1,apc2,apv1,apv2,time]=readantex('..\..\data\igs05.atx');

for i=1:length(ants)
    disp(sprintf('%-20s: %7.4f %7.4f %7.4f  %7.4f %7.4f %7.4f',ants{i},apc1(i,:),apc2(i,:)))
end

k=find(strcmp('LEIAT504GG',ants));
%k=find(strcmp('ASH701945C_MSCIS',ants));
az=0:5:360;
ze=0:5:90;

for i=1:length(az)
    msg=sprintf('%3.0f: ',az(i));
    for j=1:length(ze)
        msg=[msg,sprintf(' %6.2f',apv1(k,i,j)*1E3)];
    end
    disp(msg);
end
disp('');
for i=1:length(az)
    msg=sprintf('%3.0f: ',az(i));
    for j=1:length(ze)
        msg=[msg,sprintf(' %6.2f',apv2(k,i,j)*1E3)];
    end
    disp(msg);
end

