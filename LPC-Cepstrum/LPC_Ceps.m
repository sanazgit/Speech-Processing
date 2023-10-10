clear;clc;

[x,Fs]=audioread('One Phone - Vaakdaar - TIMIT-MCCS0.wav');

% L=length(x);
N_fft= 1024; % Size of FFT
N=450; 	% Frame Size  ( Length )

x=x./max(abs(x)); % Normalised Signal

x1= x(1:N);
xh=x1.*hanning(N);
a_spec= fft(xh,N_fft);
H1= a_spec/max(a_spec);

f= (0:N_fft/2)*Fs/N_fft;
plot(f,log(eps+abs(H1(1:N_fft/2+1)))); hold on

% lpc Analysis %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[a,e]=lpc(x,24);
b0= sqrt(e);
b=(b0);
Hz= freqz(b,a,4000);

H2= abs(Hz)/max(abs(Hz));

plot(log(H2),'r'); hold on

% g= filter(a,b,x);
% g_spec= fft(g,N_fft);
% Hg= g_spec/max(g_spec);
% plot(f,log(eps+abs(Hg(1:N_fft/2+1))),'y')

% Cep Analysis %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
cx = real(ifft(eps+log10(a_spec),N_fft));
L=12;
win= hanning(2*L+1)';
win= win/max(win);
lifter= [ones(1,L+1) zeros(1,N_fft-(2*L+1)) ones(1,L)]';
% lifter= [1 win(L+2:2*L+1) zeros(1,N_fft-(2*L+1)) win(1:L)]';
clifter=  cx .*lifter;

Xk= fft(clifter,N_fft);
Xk= Xk-max(real(Xk));
plot(f,real(Xk(1:N_fft/2+1)),'k');

grid minor
title('Absolute FFT plot')
xlabel('Frequency(HZ)')
ylabel('Log spectrum (dB)')
hold off;
legend('Signal spectrum', 'LPC Analysis', 'Inverse filter (g)','CEPS Analysis','Location', 'southwest')

