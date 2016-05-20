tic
clear all; clc;

N = 6; g =1; 
A = 0.03; 
MomValue = [ -1 0 0 1 1 2 2 3 3 4 4 5 5 6 6 7 7 8 8]'; 
LLVal = 1;
AllBras = csvread('StatesLLL.csv');

AllKets = AllBras';

LengthAllBras = length(AllBras);
LengthVect = length(AllBras(1,:));
NMatLLL(N, LengthAllBras); %Function
LMatLLL(LLVal); % Function
DeltaInteractionVariableLL(LLVal, LengthVect);%Function
InteractionTermLLL(g, LengthVect, LLVal);%Function
PotentialLLL(LengthVect, LLVal); %Function



ss=0;
aa=0; %aa is for eigenvector storage loop
%0.6:0.001:1.2 for graphing purposes 0.776 for SPDM
for Omega = 0.7:0.001:0.8
    ss=ss+1
    aa = aa+1;
    
    NAll = csvread('NMatLLLFinal.csv');
    UAll = csvread('UMatFinalLLL.csv');
    LAll = (csvread('LMatLLL.csv')).*(1-Omega);
    VNoA = csvread('VMatNoA_LLL.csv');
    VAll = VNoA.*0.5.*A;
    
    Total = NAll+UAll+LAll+VAll;
    Eig = eig(Total);
    [V, D] = eig(Total);
    TESTEigenMatrixOmega(:,aa) = Omega;
    TESTEigenMatrix1(:,aa) = V(:,1);
    OmegaIncreaseMatrix(1,ss) = Omega;
    OmegaEnergyMatrix(:,ss) = Eig;
    OmegaMatRestricted = OmegaEnergyMatrix(1:9,:);
    
end

plot(OmegaIncreaseMatrix,OmegaMatRestricted);
xlabel('Omega')
ylabel('<E>')
ylim([8.2 8.9])
xlim([0.7 0.8])
hold on; 
     
toc

