clear;clc;

[x,Fs]=audioread('S21---10sec-Male---Clean.wav');

L_total=length(x); % Total signal length
FrameShift=0.010;
FrameSize=0.010;
N= floor(FrameSize * Fs); 	% Frame Size  ( Length )
R= floor(FrameShift * Fs);	% Frame Shift ( Step )
M=floor( (L_total-N)/R + 1 ); % Number of FramesS21_NewClean

Alpha=0.7;
Beta=0.975;
Emean=0.1;
Edev= 0.25*Emean;
 
Em= zeros(1,M-1);
Threshold= zeros(1,M-1);
E_smooth= zeros(1,M-1); 
X_spec= zeros(N,M-1);

for m=1:M-1
    
   x_one_frame = x(1+(m)*R : N+(m)*R );
   win= hanning(length(x_one_frame));
   Xk= fft((x_one_frame.*win),N);
   
   Xk_db= 20*log10(eps + abs(Xk));
   X_spec(:,m)= Xk_db;
  
   thr= 1.1*(Emean + 2*Edev);
   Threshold(1,m)=thr;
   
   em= Energy(Xk,N);
   
     Em(1,m)=em;
   
   if m==1
       Es= em;
   else
       Es= Alpha*Es+(1-Alpha)*em;
   end
   
   E_smooth(1,m)= Es;
   
   if Es>= thr
       flag=1;
       Voice(:,m)= ifft(Xk);
   else
       flag=0;
   end
   
   if flag==0
       Emean= Beta * Emean + (1-Beta)*Es;
       Edev= Beta * Edev + (1-Beta)* abs(Es - Emean);
   end
   
end

figure(1)
subplot(2,1,1)
imagesc(X_spec); title('2D plot of STFT')

xlabel('Time')
ylabel('Frequency')

subplot(2,1,2)
mesh(X_spec); title('3D plot of STFT')
xlabel('Time')
ylabel('Frequency')
zlabel('Power of Frequency')

%======================= PLotting
figure(2)
plot(Em,'k'); hold on
plot(E_smooth,'m'); hold on
plot(Threshold, 'r')

legend('Energy', 'E_smooth', 'Threshold','Location', 'bestoutside')

% Voice=reshape(Voice,1,[]);
% Voice=Voice(:);
% 
% Voice= Clean(Voice);
% 
% audiowrite('S21_NewNoisy_5.wav',Voice,Fs); %Save as new sound without silent


%**************** Mean Functions
function [out]= mean(x,N)
    
    sum=0;
    
    for n=1:N
        sum= sum + x(n,1);
    end

    out= sum/N;
end

%*************** Energy Function
function [out]= Energy(x,N)

    x_mean= mean(x,N);
    
    x_new= x - x_mean;
    sum=0;
    
    for n=1:N
        sum= sum + abs(x_new(n,1));
    end

    out= sum/N;
end

% ***************** Clean Voice
function [out]= Clean(x)

    [m,n]= size(x);
    
    for i=1:m
        if x(i,1)~=0
            out(i,1)=x(i,1);
        end
    end
    
end

