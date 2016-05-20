% Changes - NEW LLVAL -   1 for LLL only, 2 for LLL and LL2
% Run once then comment out lines 17 - 23, then can alter A very quickly
% without having to re-run everything MAKE SURE DELETEMATRIX = 0
tic
clear all; clc;
DeleteMatrix = 0; % 1 = yes
N = 6; g =1; LLVal = 1;% For the deltainteractionfunction, LLL = 1, LL2incLLL = 2
A = 0.03; 
MomValue = [ -1 0 0 1 1 2 2 3 3 4 4 5 5 6 6 7 7 8 8]'; 
if LLVal == 1
AllBras = csvread('StatesLLL.csv');
elseif LLVavl ==2
    AllBras = csvread('StatesLLL_and_LL2.csv');
end
AllKets = AllBras';

LengthAllBras = length(AllBras);
LengthVect = length(AllBras(1,:));
NMatLL2(N, LengthAllBras); %Function
LMatLL2(LLVal); % Function
DeltaInteractionVariableLL(LLVal, LengthVect);%Function
InteractionTermLL2(g, LengthVect, LLVal);%Function
PotentialLL2(LengthVect, LLVal); %Function

toc

ss=0;
for Omega = 0.6:0.001:1.2
    ss=ss+1;
    
    NAll = csvread('NMatLL2Final.csv');
    UAll = csvread('UMatFinalLL2.csv');
    LAll = (csvread('LMatLL2.csv')).*(1-Omega);
    VNoA = csvread('VMatNoA_LL2.csv');
    VAll = VNoA.*A;
    
    Total = NAll+UAll+LAll+VAll;
    Eig = eig(Total);
    
    OmegaIncreaseMatrix(1,ss) = Omega;
    OmegaEnergyMatrix(:,ss) = Eig;
    OmegaMatRestricted = OmegaEnergyMatrix(1:9,:);
    
end

plot(OmegaIncreaseMatrix,OmegaMatRestricted);
xlabel('Omega')
ylabel('<E>')
ylim([8 9.2])
xlim([0.6 0.85])
hold on;

if DeleteMatrix == 1
    delete('DeltaInteractionVariableLL.csv','LMatLLL.csv','NMatLLLFinal.csv','UMatFinalLLL.csv','VMatNoA_LLL.csv')
end    
     


