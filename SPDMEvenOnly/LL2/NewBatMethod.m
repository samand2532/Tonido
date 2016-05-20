%%%% NewBat calculation, using selection to create bat states - so instead
%%%% of |6,0> after OmegaC showing as m=1, it collects these and forces
%%%% them into the |0,6> column.


%%% SPDMmk2 and Expansion Combined

clc; clear all;
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

%% Hamiltonion Omega
for Omega = 0.8285;%0.001:0.9
    Omega
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
end

% figure
% plot(OmegaIncreaseMatrix,OmegaMatRestricted);
%  xlabel('Omega')
%  ylabel('<E>')
%%SPDM
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
%% EndSPDM
%% Finds Lasgerst coeff of largest and sec largest EigVect
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
%%%% finds the 3 most populated states of the largestVect

LargestVectTemp(:,1) = LargestVect(:);
LargestVectTemp(:,2) = abs(LargestVect(:));
SecLargeVectTemp(:,1) = SecLargeVect(:);
SecLargeVectTemp(:,2) = abs(SecLargeVect(:));

Vects = [LargestVect SecLargeVect]

[C0Valmax,C0Pos] = max(LargestVectTemp(:,2));
C0Valmax = LargestVectTemp(C0Pos,1) ;
LargestVectTemp(C0Pos,1) = 0;
LargestVectTemp(C0Pos,2) = 0;

[C2Valmax,C2Pos] = max(LargestVectTemp(:,2));
C2Valmax = LargestVectTemp(C2Pos,1) ;
LargestVectTemp(C2Pos,1) = 0;
LargestVectTemp(C2Pos,2) = 0;

[C4Valmax,C4Pos] = max(LargestVectTemp(:,2));
C4Valmax = LargestVectTemp(C4Pos,1) ;
LargestVectTemp(C4Pos,1) = 0;
LargestVectTemp(C4Pos,2) = 0;

[C1Valmax,C1Pos] = max(SecLargeVectTemp(:,2));
C1Valmax = SecLargeVectTemp(C1Pos,1) ;
SecLargeVectTemp(C1Pos,1) = 0;
SecLargeVectTemp(C1Pos,2) = 0;

[C3Valmax,C3Pos] = max(SecLargeVectTemp(:,2));
C3Valmax = SecLargeVectTemp(C3Pos,1) ;
SecLargeVectTemp(C3Pos,1) = 0;
SecLargeVectTemp(C3Pos,2) = 0;

EvenFidelity = C0Valmax.^2 + C2Valmax.^2 + C4Valmax.^2
OddFidelity = C1Valmax.^2 + C3Valmax.^2 

C0Val = C0Valmax;
C2Val = C2Valmax;
C4Val = C4Valmax;
C1Val = C1Valmax;
C3Val = C3Valmax;

