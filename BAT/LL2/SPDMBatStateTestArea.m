%%%% SPDMmk2 and Expansion Combined
clc; clear all; 
tic
A = 0.03;
g = 1;
kk=0;
N = 6;
Basis = csvread('EvenBasisLL2.csv');
AllBras = Basis;
AllKets = AllBras';

NMat = csvread('NMat.csv');
absLMat = csvread('absLMat.csv');
littlenMat = csvread('littlenMat.csv');

UMat = csvread('MyUMat.csv');
LMat = csvread('LMat.csv');
VMatnoA = csvread('VMatnoA.csv');

ss=0;

 for Omega = 0.82856
    
    ss=ss+1;
    
    NAll = NMat;
    LAll = LMat;
    absL = absLMat;
    
    UAll = g.*UMat;
    VAll = VMatnoA.*A;
   
    littlenAll = 2.*littlenMat;
    Total = littlenAll+(absLMat-(Omega*LAll))+(NAll)+UAll+VAll;
    Eig = eig(round(Total,15));
    
    OmegaIncreaseMatrix(1,ss) = Omega;
    OmegaEnergyMatrix(:,ss) = Eig;
    OmegaMatRestricted = OmegaEnergyMatrix(1:9,:);
    
    [VHam,DHam] = eig(Total);


% figure
% plot(OmegaIncreaseMatrix,OmegaMatRestricted);
%  xlabel('Omega')
%  ylabel('<E>')
Dordered = diag(DHam);

[HamMinVal,HamMinPos] = min(Dordered(:));


[x,y] = size(Basis);
SPDMmat = zeros(y);
creaVect = zeros(y,1);
annVect = zeros(y,1);

OrigCoeff  = VHam(:,HamMinPos);

% Norm = 1/sqrt(192);
% OrigEigValNormalised = OrigEigVal.*Norm;

for Kets = 1:x
    for k = 1:y
        for l = 1:y
            C1 =OrigCoeff(Kets,1);
            
            InitialKet = AllKets(:,Kets);
    creaVect = creaVect.*0;
    annVect = annVect.*0;
    creaVect(k,1) = 1;
    annVect(l,1) = -1;
    
    annConst = sqrt(InitialKet(l,1));
    annTrans = InitialKet + annVect;
    creaConst = sqrt(annTrans(k,1)+1);
    FinalKet = annTrans + creaVect;
    
    for Bras = 1:x
        Bra = AllBras(Bras,:);
        C2 =OrigCoeff(Bras,1);
        if Bra == FinalKet'
            SPDMmat(k,l) = SPDMmat(k,l) + (creaConst*annConst*C1*C2);
        end
    end
        end
    end
end


[eigVectSPDM,eigValSPDM] = eig(SPDMmat);

%%%%% finds 2 largest (+ve) eigenstates
%%% by diaganolising the eigValSPDM mat, finding the 
%%% corresponding row, then overwriting it with 0 then
%%% finding the next max.

diageigValSPDM = diag(eigValSPDM);
[maxeigValSPDM,maxPosSPDM] = max(diageigValSPDM);
diageigValSPDM(maxPosSPDM,1) = 0;
[secmaxeigVectSPDM,secmaxPosSPDM] = max(diageigValSPDM);

[maxeigValSPDM,maxPosSPDM];
[secmaxeigVectSPDM,secmaxPosSPDM];

LargestVect = eigVectSPDM(:,maxPosSPDM);
SecLargeVect = eigVectSPDM(:,secmaxPosSPDM);

C2 = LargestVect(2,1)
C4 = LargestVect(4,1)
C13= LargestVect(13,1)
C1 = SecLargeVect(1,1)
C3 = SecLargeVect(3,1)
Psi1before = (C2^2)+(C4^2)+(C13^2)
Psi2before = (C1^2)+(C3^2)


C2a = SecLargeVect(2,1)
C4a = SecLargeVect(4,1)
C13a= SecLargeVect(13,1)
C1a = LargestVect(1,1)
C3a = LargestVect(3,1)
Psi1after = (C2a^2)+(C4a^2)+(C13a^2)
Psi2after = (C1a^2)+(C3a^2)
kk = kk + 1
Table(1,kk) = Omega;
Table(2,kk) = Psi1before;
Table(3,kk) = Psi2before;
Table(4,kk) = Psi1after;
Table(5,kk) = Psi2after;




 end
toc