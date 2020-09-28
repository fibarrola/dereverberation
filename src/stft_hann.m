function spc=stft_hann(x,nfft,win,over)

% Returns short time fourier transform of signal x
% spc=stft(x,nfft,win,over)
%
% x      signal - must be a row vector
% nfft   length to use for the fourier transform (can be used for zero-padding).
% win    length of hamming window to use. Defines the frame length, usually equal to nfft
% over   number of overlaping samples between two consecutive frames.
%
% The size of the transform is nfft by fix((nx-over)/(length(win)-over))
% Author: Leandro Di Persia
% Copyright: Doshisha University

[m,n]=size(x);
%check if it is row vector
nwin=win;
win=hann(nwin);
% win=hamming(nwin);
%win=boxcar(nwin);
if n==1,
    x=x.';
end
[m,n]=size(x);
%calculate number of columns in the matrix
ncol=fix((n-over)/(nwin-over));
mlen=ncol*(nwin-over)+over;
if n < mlen,
    x(mlen)=0;
else
    if n > mlen
       mlen=mlen+nwin-over;
       ncol=ncol+1;
       x(mlen)=0;
    end
end
xmat=zeros(nwin,ncol);
for f=1:ncol,
    init=(f-1)*(nwin-over)+1;
    a=x(init:init+nwin-1);
    xmat(:,f)= a.';
end
for k=1:ncol,
    xmat(:,k)=xmat(:,k).*win;
end
spc=fft(xmat,nfft);

     