clear;clc;
 
[x,Fs]=audioread('a_vac.wav');
 
L_total=length(x);
N_fft= 2048; % Size of FFT
N=450;  % Frame Size  ( Length )
 
x=x./max(abs(x)); % Normalised Signal
 
x1= x(1:N);
xh=x1.*hanning(N);
a_spec= fft(xh,N_fft);
 
f= (0:N_fft/2)*Fs/N_fft;
 
% plot(log10(abs(a_spec)))
plot(f, 20*log10(abs(a_spec(1:N_fft/2+1))));
 
grid minor
title('Absolute FFT plot')
xlabel('Frequency(HZ)')
ylabel('Log spectrum (dB)')
