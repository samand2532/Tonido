clc; clear all; close all;
tic 
%%%% Recaluculate Terms?
ReCalIntTerm = 0; 
ReCalAnisoTerm = 0; %%%% 0 = no, 1 = yes

N = 6;
g = 1;
A = 0.03;
LLVal=1;
KPos = [ 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19; %Matlab position
         0 0 1 0 1 0 1 0 1 0  1  0  1  0  1  0  1  0  1 ; %n value
        -1 0 0 1 1 2 2 3 3 4  4  5  5  6  6  7  7  8  8]; %Mt value
    
AllBras = csvread('NewLL2.csv');
AllKets = AllBras';
LengthBasis = length(AllBras);

NMat = eye(LengthBasis).*N;

LMat = zeros(LengthBasis);

for LBra = 1:LengthBasis
    InitialBraL = AllBras(LBra,:);
    for LKet = 1:LengthBasis
        InitialKetL = AllKets(:,LKet);
        
        if InitialBraL == InitialKetL'
            LMat(LBra,LKet) = dot(InitialBraL,KPos(3,:));
        end
    end
end
%%%% InteractionTerm area
if ReCalIntTerm == 0
   UMatNoConst = csvread('UMatNoConst.csv');
elseif ReCalIntTerm == 1
   IntFunc();
end
%%%%%%%%%%%%%%%%%%%%%%%%

%%%% Anisotropic Term area
if ReCalAnisoTerm == 0
    AnisoTermNoA = csvread('AnisotermNoConst.csv');
elseif ReCalAnisoTerm == 1
    AnisoFunc();
end
%%%%%%%%%%%%%

UMatNoConst = csvread('UMatNoConst.csv');
UMat = UMatNoConst.* (g/(4*pi));
AnisoTermNoA = csvread('AnisotermNoConst.csv');
AnisoMat = A.*AnisoTermNoA;



ss=0;
for Omega = 0.6:0.01:0.9
    ss=ss+1;
    
    NAll = NMat;
    UAll = UMat;
    LAll = LMat.*(1-Omega);
    VAll = 0.5*AnisoMat;
    
    Total = NAll+UAll+LAll+VAll;
    Eig = eig(Total);
    [V,D] = eig(Total);
    
    OmegaIncreaseMatrix(1,ss) = Omega;
    OmegaEnergyMatrix(:,ss) = Eig;
    OmegaMatRestricted = OmegaEnergyMatrix(1:9,:);
    
end

plot(OmegaIncreaseMatrix,OmegaEnergyMatrix);
xlabel('Omega')
ylabel('<E>')
ylim([7.7 8.5])
xlim([0.74 0.86])
hold off;








    toc
 