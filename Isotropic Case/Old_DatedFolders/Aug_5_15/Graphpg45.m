%%%%%%% Bad code but able to re-create fig 3.1, need to manually 
%adjust L, ni and Negmi

clc; clear all; 

N = 6;
ni=0;
NegMi = 0;
L=0;
i=1;
    
    for Omega = 0:0.001:1
        i=i+1;
        Storage(i,1) = Omega;
        E = L.*(1-Omega) + 2.*(NegMi+ni);
        Storage(i,2) = E;
    end
    
plot(Storage(:,1),Storage(:,2))
xlim([0.001 1])
hold on;