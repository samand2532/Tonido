clc; clear all

tic
A = 0.03;
g = 1;

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
format Long
asd=0;
%% Hamiltonion Omega
Omega = 0.83
    %0.818567080538138;
    
    
    %%%Omega = 0.825805 for roughly bat
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

Dordered = diag(DHam);
[HamMinVal,HamMinPos] = min(Dordered(:));
[x,y] = size(Basis);
SPDMmat = zeros(y);
creaVect = zeros(y,1);
annVect = zeros(y,1);
OrigCoeff  = VHam(:,HamMinPos);
Wavefunction = [OrigCoeff Basis];

for Kets = 1:x
    for k = 1:y
        for l = 1:y
            CL1 =OrigCoeff(Kets,1);
            
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
                CL2 =OrigCoeff(Bras,1);
                if Bra == FinalKet'
                    SPDMmat(k,l) = SPDMmat(k,l) + (creaConst*annConst*CL1*CL2);
                end
            end
        end
    end
end
[eigVectSPDM,eigValSPDM] = eig(SPDMmat);

diageigValSPDM = diag(eigValSPDM);
[maxeigValSPDM,maxPosSPDM] = max(diageigValSPDM);
diageigValSPDM(maxPosSPDM,1) = 0;
[secmaxeigVectSPDM,secmaxPosSPDM] = max(diageigValSPDM);

[maxeigValSPDM,maxPosSPDM];
[secmaxeigVectSPDM,secmaxPosSPDM];
LargestVect = eigVectSPDM(:,maxPosSPDM);
SecLargeVect = eigVectSPDM(:,secmaxPosSPDM);
Vects = [LargestVect SecLargeVect]

%%% finds max val and position for Gs and 1st ES

GSTemp(:,1) = abs(Vects(:,1));
GSTemp(:,2) = Vects(:,1);
Ex1(:,1) = abs(Vects(:,2));
Ex1(:,2) = Vects(:,2);

[GSTVal,GSTPos] = max(GSTemp(:,1));
GSVal = GSTemp(GSTPos,2);
GSPos = GSTPos;

[Ex1TVal,Ex1TPos] = max(Ex1(:,1));
Ex1Val = Ex1(Ex1TPos,2);
Ex1Pos = Ex1TPos;

[GSVal GSPos]
[Ex1Val Ex1Pos]
































toc