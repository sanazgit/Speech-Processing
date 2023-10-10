clear;close all;clc;

[x,Fs]=audioread('One Phone - Vaakdaar - TIMIT-MCCS0.wav');

N=256; 	% Frame Size 

seg= x(1:N);
wseg=seg.*hanning(N);

P = 16;
[a,G] = lpc(wseg,P);

c0= log(G);
a=a(2:end);

for n=1:N-1
    
    if n==1
        Cn(1,1)=a(1);
    end
    
    if n>1 && n<=P
        sum=0;
        for k=1:n-1
            sum= sum + ((k/n)*Cn(k)*a(n-k));   
        end
        Cn(1,n)= -(a(n) + sum);
    end    
    
     
    if n>P
        sum=0;
        for k=n-P:n-1
            sum= sum + ((k/n)*Cn(k)*a(n-k));   
        end
        
        Cn(1,n)= -sum;
        
    end

end

Cn=[c0 Cn];

subplot(2,1,1)
plot(a); title('LPC Coefficients')

subplot(2,1,2)
plot(Cn(1:100)); title('LPCC Coefficients')

