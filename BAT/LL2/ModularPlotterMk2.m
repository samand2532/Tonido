%%%% Modular approach to graphing
%%%% THIS IS ADDED AS A GIT TEST

clc; clear all;
tic
A = 0.03;
g = 1;
%AllBras = csvread('LuisBasisLL2.csv');
AllBras = csvread('EvenBasisLL2.csv');

NMat = csvread('NMat.csv');
absLMat = csvread('absLMat.csv');
littlenMat = csvread('littlenMat.csv');
%UMat = csvread('LuisU.csv');
UMat = csvread('MyUMat.csv');
LMat = csvread('LMat.csv');
VMatnoA = csvread('VMatnoA.csv');


ss=0;
for Omega = 0.83%0.7:0.0001:0.95
    
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

%plot(OmegaIncreaseMatrix,OmegaMatRestricted);
% xlabel('Omega')
% ylabel('<E>')
% 
%   ylim([7.7 8.5])
%   xlim([0.74 0.95])
 

toc