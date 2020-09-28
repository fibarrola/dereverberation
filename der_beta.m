% A blind, single channel dereverberation method in the time-frequency 
% domain that makes use of variable beta-divergence as a cost function.
% 
% The algorithm needs a .wav file corresponding to the reverberant signal
% (isrev = 1), or a clean signal that it will make reverberant if isrev = 0.
% 
% Parameters
% wsize: window size (in samples) for STFT.
% overl: window overlapping (in samples) for STFT.
% M: Time frames for RIR matrix H.
% J: Number of dictionay elements.
% beta1: parameter for beta-divergence in Stage 1 (dictionary learning).
% beta2: parameter for beta-divergence in Stage 2 (representation).
% ldaU: penalization parameter for building U in Stage 2.
% ldaH: penalization parameter for building H in Stage 2.

% Recommended values for 16[kHz] signals are: wsize=512, overl=256, M=20, 
% J=64, beta1=0.75, beta2=2, ldaU=0.1, ldaH = 0.3. For artifficial 
% (noiseless) dereverberation, it is recommended to set ldaU=1e-4;

function [der_sig,rev_sig,specs,rep] = der_beta(sig_name,isrev,wsize,overl,M,J,beta1,beta2,ldaU,ldaH)

% Preliminaries -----------------------------------------------------------
addpath(genpath(pwd));
load colores

% Default settings --------------------------------------------------------
if nargin == 2
    p.wsize = 512;
    p.overl = 256;
    p.M = 20;
    p.J = 64;
    p.beta1 = 0.75;
    p.beta2 = 2;
    p.ldaU = 0.0001+0.1*isrev;
    p.ldaH = 0.3;
else
    p.wsize = wsize;
    p.overl = overl;
    p.M = M;
    p.J = J;
    p.beta1 = beta1;
    p.beta2 = beta2;
    p.ldaU = ldaU;
    p.ldaH = ldaH;
end
p.it1 = 250;
p.it2 = 20;
p.delta = 1e-7;


% Building data -----------------------------------------------------------
if isrev,
    [s_rev,fs] = wavread(sig_name);
else
    [s_cln,fs] = wavread(sig_name);
    s_rir = wavread('rir600ms16khz.wav');
    s_rev = filter(s_rir,1,s_cln);
    s_rev = s_rev/max(abs(s_rev));
end

% Preprocessing -----------------------------------------------------------
b = fir1(5000,60/fs,'high'); % only if records
s_rev = conv(s_rev,b);
s_rev = s_rev(2501:end-2500);

% Making spectrograms -----------------------------------------------------
[Y_rev,N,K] = makespec_hann(s_rev,fs,p.wsize,p.wsize,p.overl);
Y = Y_rev.*conj(Y_rev);

% Variable beta-divergence NMF --------------------------------------------
[W,U,H,S] = fact_beta(Y,p);

% Reconstruction ----------------------------------------------------------
F = sqrt(S).*exp(1i*atan2(imag(Y_rev),real(Y_rev)));
aux   = flipud(conj(F(2:K-1,:)));
s_der = istft_hann([F;aux],p.wsize,p.wsize,p.overl)';
s_der = s_der/max(abs(s_der));

% Plotting and saving spectrograms ----------------------------------------
figure,
subplot(2,1,1), imagesc((0:N-1)/2^7,(0:256)*8/256,log(Y+1e-2)), axis xy, xlabel('time [s]'), ylabel('frequency [kHz]'),title('Reverberant spectrogram Y');
subplot(2,1,2), imagesc((0:N-1)/2^7,(0:256)*8/256,log(S+1e-2)), axis xy, xlabel('time [s]'), ylabel('frequency [kHz]'),title('Restored spectrogram S');
colormap(gr2bl)
specs = ['out' filesep 'spectrograms.jpg'];
fileout2 = [pwd filesep 'out' filesep 'spectrograms.jpg'];
print('-djpeg90',fileout2);
pause(0.1);

% Plotting and saving representation elements -----------------------------
figure,
[W,U] = sortvecs(W,U);
subplot(3,4,[2,3]), imagesc(log(U+1e-2)); axis xy, title(sprintf('Coefficient matrix U'));
subplot(3,4,[5,9]), imagesc(1:p.J,(0:256)*8/256,log(W+1e-2)); axis xy, ylabel('frequency [kHz]'), title('Dictionary W');
subplot(3,4,[6,7,10,11]),imagesc((0:N-1)/2^7,(0:256)*8/256,log(S+1e-2)), axis xy, xlabel('time [s]'), title('Restored spectrogram S');
subplot(3,4,[8,12]),imagesc((0:p.M-1)/2^7,(0:256)*8/256,log(H+1e-2)), axis xy, xlabel('time [s]'), title('RIR H');
colormap(gr2bl)
rep = ['out' filesep 'representation.jpg'];
fileout3 = [pwd filesep 'out' filesep 'representation.jpg'];
print('-djpeg90',fileout3);

% Saving ------------------------------------------------------------------
der_sig = ['out' filesep 's_der.wav'];
fileout1 = [pwd filesep 'out' filesep 's_der.wav'];
wavwrite(s_der,fs,fileout1);

rev_sig = ['out' filesep 's_rev.wav'];
fileout1 = [pwd filesep 'out' filesep 's_rev.wav'];
wavwrite(s_rev,fs,fileout1);

end

