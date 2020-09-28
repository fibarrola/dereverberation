% [Y,T,X] = makespec(signal,fs,wsize,wsize2,overl) makes a spectrogram Y of
% size T by X for the specified the signal, with sampling frequency fs.
% wsize determines the window size and overl the overlapping.
%
function [Y,T,X] = makespec_hann(signal,fs,wsize,wsize2,overl) 
tt = 0:1/fs:(length(signal)-1)/fs;    % Number of signal dots
Y  = stft_hann(signal,wsize,wsize2,overl); % Short-time Fourier transform
Y  = Y(1:size(Y,1)/2+1,:);            % Taking A matrix's upper half
[X,T] = size(Y);
end