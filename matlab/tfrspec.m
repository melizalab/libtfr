function SP = tfrspec(S, N, step, Np, K, tm, flock, tlock)
% TFRSPEC  Compute multitaper time-frequency reassignment spectrogram
% 
% Usage:
%    SP = TFRSPEC(S, N, step, Np, K, tm, flock, tlock)
%
%         S - input signal (short integers, floats, doubles are supported)
%         N - number of points in FFT window
%         step - number of points to shift for each spectrogram frame
%         Np - number of points in tapers
%         K  - number of tapers (default 6)
%         tm - time support of tapers (default 6.0)
%         flock - frequency locking (in Hz; default 0.01)
%         tlock - time locking (in frames; default 5)
%
%         SP - N/2+1 by length(S)/step power spectrogram (linear
%         scale)
%