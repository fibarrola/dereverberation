function x=istft_hann(spc,nfft,win,ovr)

% x=istft(spc,nfft,ovr,win)
% Takes inverse short time fourier transform of spc.
%
% spc    time-frequency spectrum to reconstruct
% nfft   number of samples to use with ifft.
% ovr    number of samples to overlap between adyacent segments of signal.
% win    length of hamming window used in transform

[n,m]=size(spc);
nwin=win;
win=hann(nwin);
% win=hamming(nwin);
%win=boxcar(nwin);

nx = m*(nwin-ovr)+ovr;
x=zeros(1,nx);

xmat=ifft(spc,nfft);
ven=zeros(1,nx);
for k=1:m,
    init =(k-1)*(nwin-ovr)+1;
    a=xmat(1:nwin,k);
    b=length(a);
    ven(init:init+nwin-1)=ven(init:init+nwin-1)+win';
    x(init:init+nwin-1)=x(init:init+nwin-1)+a.';
end;
ven(1:ovr)=ven(ovr+1);
ven(nx-ovr+1:end)=ven(ovr+1);
x=real(x./ven);
%x=real(x);