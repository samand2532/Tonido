clear all;
hbar = 6.63e-34/(2*pi);
g=1;
n=0;
Wp=0.75;

for Omega = 0.0:0.001:0.85;
    n=n+1;
    
    NAll = eye(39).*6;
    V_003 = csvread('V003AllStates.csv');
    UAll = csvread('UmatAll.csv');
    LAllread = csvread('LMatAllExcludingOmega.csv');
    LAll = LAllread.*(1-Omega);
    
    
    Total = NAll+(0.5*V_003)+UAll+LAll; % With V=0.03
    Eig = eig(Total)
    
 %    Total = NAll+UAll+LAll; % With V=0.00
  %   Eig = eig(Total)

    
    OmegaIncreaseMatrix(1,n) = Omega;
    OmegaEnergyMatrix(:,n) = Eig;
   
    
    
end

plot(OmegaIncreaseMatrix,OmegaEnergyMatrix);
xlabel('Omega')
ylabel('<E>')
ylim([8.2 8.9])
xlim([0.7 0.80])
hold on;
