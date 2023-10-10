clear;close all;clc;

[x,Fs]=audioread('mic_F01_si474-11kHz---for-Pitch.wav');

timeAxis=(1:length(x))/Fs; % convert sample number to time 

figure(1)
plot(timeAxis,x); 
title('Waveform');
xlabel('Time(seconds)');

x= x/max(x);
L_total=length(x); % Total signal length

FrameShift=0.01;
FrameSize=0.032;
N= floor(FrameSize * Fs);  	% Frame Size  ( Length )
R= floor(FrameShift * Fs);	% Frame Shift ( Step )
M=floor( (L_total-N)/R + 1 ); % Number of FramesS21_NewClean

F0_m= zeros(1,M-1);

kmin= floor(Fs/600);      
kmax=floor(Fs/50);
mf = dsp.MedianFilter(2);

for m=1:M-1
    
   x_frame = x(1+(m)*R : N+(m)*R);
   x_frame_normal= x_frame - mean(x_frame);
   win= hanning(length(x_frame_normal));
   seg= x_frame_normal.*win;
   
   Rm= xcorr(seg);
   R_norm= Rm /(mean(abs(seg)));
   Nc= length(R_norm);
   R_norm= R_norm(ceil(Nc/2):end);
   
   [Am,index_L0]= max(R_norm(kmin:kmax));
   
   L0= index_L0+kmin-1;
   F0= Fs/L0;

   if (Am < 0.5)
       
       F0_m(1,m)= 0;
   else
       
       F0_m(1,m)= F0;
   end  
end

y = mf(F0_m); 

figure(2)
plot(F0_m); hold on
plot(y,'g')


