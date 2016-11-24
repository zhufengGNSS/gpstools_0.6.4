function y=bpfilt(time,x,lowf,highf)
%-------------------------------------------------------------------------------
% [system] : GpsTools
% [module] : bandpass filter
% [func]   : bandpass filter
% [argin]  : time    = time(sec)
%            x       = values
%            lowf    = low-cut frequency (hz) (0:dc)
%            highf   = high-cut frequency (hz)
% [argout] : y       = band-pass filtered values
% [note]   :
% [version]: $Revision: 3 $ $Date: 06/07/08 14:10 $
%            Copyright(c) 2004-2006 by T.Takasu, all rights reserved
% [history]: 06/05/04  0.1  new
%-------------------------------------------------------------------------------
y=repmat(nan,size(x));
i=find(all(~isnan(x),2)); if length(i)<2, return, end
x=interp1(time(i),x(i,:),time(:));
i=find(all(~isnan(x),2));
f=(0:length(i)/2)/(time(i(end))-time(i(1)));
j=find(lowf<=f&f<=highf); k=length(i)-j+2; k(k>length(i))=[];
filt=zeros(length(i),size(x,2)); filt([j,k],:)=1;
if 0<lowf, x(i,:)=detrend(x(i,:)); end
y(i,:)=real(ifft(filt.*fft(x(i,:))));
