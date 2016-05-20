clc; clear all; close all;
%%%%_____________
%%%_Variables
Nt = 6;
g = 1;

%%%%%%________________
tic
AllBras = csvread('StatesLL2.csv');
AllKets = AllBras';

KPos = [ 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19;%%K number
    -1 0 0 1 1 2 2 3 3 4  4  5  5  6  6  7  7  8  8;%%%Mt
    0 0 1 0 1 0 1 0 1 0  1  0  1  0  1  0  1  0  1];%%% n value


NMat = eye(length(AllBras)).*Nt;
LMat = eye(length(AllBras));

for Bras = 1:length(AllBras)
    InitialBra = AllBras(Bras,:);
    for Kets = 1:length(AllKets)
        InitialKet = AllKets(:,Kets);
        if InitialBra == InitialKet';
            LMat(Bras,Kets) = dot(InitialBra,KPos(2,:));
        end
    end
end

UmatNoConst = csvread('UMatNoConst.csv');
UMat = UmatNoConst.*(g/(4*pi));




    ss=0;
for Omega = 0.6:0.1:1
    ss=ss+1;
    
    NAll = NMat;
    UAll = UMat;
    LAll = LMat.*(1-Omega);
    %VNoA = csvread('VMatNoA_LL2.csv');
    %VAll = VNoA.*A;
    
    Total = NAll+UAll+LAll;%+VAll;
    Eig = eig(Total);
    
    OmegaIncreaseMatrix(1,ss) = Omega;
    OmegaEnergyMatrix(:,ss) = Eig;
    OmegaMatRestricted = OmegaEnergyMatrix(1:9,:);
    
end

plot(OmegaIncreaseMatrix,OmegaEnergyMatrix);
xlabel('Omega')
ylabel('<E>')
ylim([8 9.2])
xlim([0.6 0.85])
hold on;















toc

            
            
            
            
