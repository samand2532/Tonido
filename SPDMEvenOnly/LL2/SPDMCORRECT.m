clc; clear all; 
tic
A = 0.03;
g = 0.6;

Nt = 6;
Basis = csvread('EvenBasisLL2.csv');
%Basis = csvread('LuisBasisLL2.csv');
AllBras = Basis;
AllKets = AllBras';

NMat = csvread('NMat.csv');
absLMat = csvread('absLMat.csv');
littlenMat = csvread('littlenMat.csv');

UMat = csvread('MyUMat.csv');
LMat = csvread('LMat.csv');
VMatnoA = csvread('VMatnoA.csv');
%VMatnoA = csvread('LuisV192.csv');


ss=0;
asd=0;
    Omega = 0.8
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
    
    [V,D] = eig(Total);


Dordered = diag(D);

[HamMinVal,HamMinPos] = min(Dordered(:));


[x,y] = size(Basis);
SPDMmat = zeros(y);
creaVect = zeros(y,1);
annVect = zeros(y,1);

 OrigCoeff  = V(:,HamMinPos);
%OrigEigVal = D;
%Norm = 1/sqrt(192);
%OrigEigValNormalised = OrigEigVal.*Norm;
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
format short

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

Vects = [LargestVect SecLargeVect];



toc