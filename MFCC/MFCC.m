%% %Main Code
clear;close all;clc;


[x,Fs]=audioread('S219_Female_FarsDat_OneSentence-8kHz-without Silence.wav');

timeAxis=(1:length(x))/Fs; % convert sample number to time 

figure(1)
plot(timeAxis,x); 
title('Waveform');
xlabel('Time(seconds)');

L_total=length(x); % Total signal length

FrameShift=0.01;
FrameSize=0.032;
N= floor(FrameSize * Fs);  	% Frame Size  ( Length )
R= floor(FrameShift * Fs);	% Frame Shift ( Step )
M= floor((L_total-N)/R + 1 ); % Number of FramesS21_NewClean
Nf= 20;
Nc= 16; 
alpha= 0.975;
Fmin= 100;
Fmax= 3900;

ccs= zeros(M-1,Nc);

for m=1:M-1
    
   x_frame = x(1+(m)*R : N+(m)*R);
   xp= pre_emphasis(alpha, x_frame);
   x_frame_normal= xp - mean(xp);
   win= hanning(length(x_frame_normal));
   seg= x_frame_normal.*win;
   ccs(m,:)= mfcc_model(seg,Nf,Nc,Fmin,Fmax,Fs);
end

% save 2 rows of mfcc table as features for each sound 

% female_mfcc= ccs(1:2,:);
% save('female_mfcc','female_mfcc');

%% Plot Examples

load ('male_mfcc.mat');
load('female_mfcc.mat');
plot(male_mfcc(1,:), male_mfcc(2,:), 'k+','LineWidth', 2, 'MarkerSize', 7); hold on
plot(female_mfcc(1,:), female_mfcc(2,:), 'ko', 'MarkerFaceColor','y','MarkerSize', 7);

xlabel('X1')
ylabel('X2')
title('MFCC Plot')
legend('Male mfcc','Female mfcc', 'Location', 'northwest')


%% % ==================  MFCC Extraction ===================

function [mel]= f2mel(hz)
    mel= 2595 * log10(1+(hz/700));
end    

function [hz]= mel2f(mel)
    hz= 700 * (10.^(mel/2595)-1);
end

function band= spread_mel(hz_points, hz_c, hz_size, hz_max)

    % hz_points is an array spaced in Hz
    % hz_c is the current index
    
    band= zeros(1, hz_size);
    hz1= hz_points(max(1,hz_c-1));                      %start
    hz2= hz_points(hz_c);                               %middle
    hz3= hz_points(min(length(hz_points),hz_c+1));      %end
    
    %---
    
    for hi=1:hz_size
        hz= hi*hz_max/hz_size;
        
        if hz>hz3
            band(hi)=0;
        elseif hz >= hz2
            band(hi)= (hz3-hz)/(hz3-hz2);
        elseif hz >= hz1
            band(hi)= (hz-hz1)/(hz2-hz1);
        else
            band(hi)=0;
        end 
    end
end


function y= pre_emphasis(alpha, x)
     B = [1 alpha];
     y = filter(B,1,x);
end

function cc= mfcc_model(seg,N,M,Fmin,Fmax,Fs)

     m_low=f2mel(Fmin); % mel span lower limit
     
     if Fmax==Fs
         m_top=f2mel(Fmax/2); % mel span upper limit
     else
         m_top=f2mel(Fmax);   
     end
     
       
     mdiv=(m_top - m_low)/(N-1); % mel resolition
 
     % define an array of center frequencies
     xm= m_low:mdiv:m_top;
 
     % convert this to Hz frequencies
     xf=mel2f(xm);
 
     % quantise to the fft resolotion
     k= floor(length(seg)* (xf/Fmax));
     % take the fft o speech ...
     S= abs(fft(seg));
     

     % compute the mel filterbank
     x1= zeros(1,N);
     for xi=1:N
         band= spread_mel(xf, xi, length(S), Fmax);
         x1(xi)= sum(S .* band');
     end
     
     x= log(x1);
     
   
     cc= zeros(1, M);
     for xc=1:M
         cc(xc)= sum(x.*cos(pi*xc*((1:N) - 0.5) / N));
     end

end

