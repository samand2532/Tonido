%%%% Modular approach to graphing


clc; clear all;
tic
A = 0.03;
g = 0.2;
%AllBras = csvread('LuisBasisLL2.csv');
AllBras = csvread('EvenBasisLL2.csv');

NMat = csvread('NMat.csv');
absLMat = csvread('absLMat.csv');
littlenMat = csvread('littlenMat.csv');

UMat = csvread('MyUMat.csv');
LMat = csvread('LMat.csv');
VMatnoA = csvread('VMatnoA.csv');
%VMatnoA = csvread('LuisV192.csv');

ss=0;
for Omega = 0.7:0.001:1.1
    
    ss=ss+1;
    
    NAll = NMat;
    LAll = LMat;
    absL = absLMat;
    
    UAll = g.*UMat;
    VAll = VMatnoA.*A;
   
    littlenAll = 2.*littlenMat;
    Total = littlenAll+(absLMat-(Omega*LAll))+(NAll)+UAll+VAll;%+round(AAll,10);
    Eig = eig(round(Total,15));
    
    OmegaIncreaseMatrix(1,ss) = Omega;
    OmegaEnergyMatrix(:,ss) = Eig;
    OmegaMatRestricted = OmegaEnergyMatrix(1:9,:);
    
    [V D] = eig(Total);
end

plot(OmegaIncreaseMatrix,OmegaMatRestricted);
xlabel('Omega')
ylabel('<E>')

  %ylim([7.5 8.5])
  xlim([0.7 1.1])
  legend(['g = ',num2str(g)])
 

toc