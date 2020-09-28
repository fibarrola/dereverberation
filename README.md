# Blind Speech Dereverberation

A blind, single channel dereverberation method in the time-frequency domain that makes use of variable beta-divergence as a cost function.


# Usage

The function der_beta takes a .wav file as input corresponding to the reverberant signal (isrev = 1), or to a clean signal that it will make reverberant if isrev = 0.

Parameter description
wsize: window size (in samples) for STFT.
overl: window overlapping (in samples) for STFT.
M: Time frames for RIR matrix H.
J: Number of dictionay elements.
beta1: parameter for beta-divergence in Stage 1 (dictionary learning).
beta2: parameter for beta-divergence in Stage 2 (representation).
ldaU: penalization parameter for building U in Stage 2.
ldaH: penalization parameter for building H in Stage 2.

Recommended values for 16[kHz] signals are: wsize=512, overl=256, M=20, J=64, beta1=0.75, beta2=2, ldaU=0.1, ldaH = 0.3. For artifficial (noiseless) dereverberation, it is recommended to set ldaU=1e-4;


# References

If useful, please cite

F. J. Ibarrola, R. D. Spies and L. E. D. Persia, "Switching Divergences for Spectral Learning in Blind Speech Dereverberation," in IEEE/ACM Transactions on Audio, Speech, and Language Processing, vol. 27, no. 5, pp. 881-891, May 2019, doi: 10.1109/TASLP.2019.2901643.
