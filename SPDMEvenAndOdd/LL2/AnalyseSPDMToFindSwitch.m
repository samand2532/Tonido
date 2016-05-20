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
Rec=0;
%% Hamiltonion Omega
for Omega = 0.828567080538137:0.00000000000000001:0.828567080538138;%0.001:0.9
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

Vects = [LargestVect SecLargeVect];
Rec = Rec + 1;
RecMat(Rec,1) = Omega;
RecMat(Rec,2) = LargestVect(2,1);
RecMat(Rec,3) = LargestVect(3,1);
if round(LargestVect(2,1),4) == 0 
    Omega
end
end
toc