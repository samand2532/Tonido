clear all; clc;
L = csvread('LMat.csv');
U= csvread('UMatFinal.csv');
V = csvread('VMatFinalincA.csv');
N = eye(39).*6;
n=0;

for Omega = 0.6:0.001:0.85
    n=n+1;
    
    LOmega = L.*(1-Omega);
    Total = N+U+(LOmega)+(V);
    Eig = eig(Total);
    
    OmegaIncreaseMatrix(1,n) = Omega;
    OmegaEnergyMatrix(:,n) = Eig;
    %%% Following line selects only the lowest 8 eigenstates
    OmegaMatRestricted = OmegaEnergyMatrix(1:9,:);
    
    
end

plot(OmegaIncreaseMatrix,OmegaMatRestricted);
xlabel('Omega')
ylabel('<E>')
ylim([8 9.2])
xlim([0.6 0.85])
hold off;