clear;close all;clc;

[x,Fs]=audioread('One Phone - Vaakdaar - TIMIT-MCCS0.wav');

N=256; 	% Frame Size 

seg= x(1:N);
wseg=seg.*hanning(N);

P = 16;
[a,G] = lpc_lv(wseg,P);    
    
% ======================= Levinson_Durbin =============================

function [a, G]= lpc_lv(sig,p)
    N= length(sig);
    r= xcorr(sig);
    R_xx=r(N:end);
    [a,G]= levinson(R_xx,p);
end

function [a_lpc, G]= levinson(R,P)

    E0= R(1); 

    for i = 1:P
        
       if i==1
           k= R(i+1)/E0;
           
       else
           sum=0;
           for j=1:i-1
               sum = sum + (a(i-1,j)*R(i-j+1));
           end
           k= (R(i+1)- sum)/E;
       end   
       
       a(i,i)= k;
       
       if i==1
           E= (1- k*k) * E0;
       else    
           for j=1:i-1
               a(i,j)= a(i-1,j)- k * a(i-1,i-j); 
           end
           E= (1- k*k) * E;
       end
    
    end
  
    for i=1:P
        a_lpc(1,i)= a(P,i); 
    end

    res=0;
    for k=1:P
        res = res+ a_lpc(k)*R(k+1);
    end
    
    G= (R(1)- res)/length(R);
    
end

