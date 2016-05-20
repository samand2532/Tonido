%%% SPDM using even and odd basis LL2

clc; clear all; 
tic
A = 0.03;
g = 1;

Nt = 6;
%Basis = csvread('EvenBasisLL2.csv');
Basis = csvread('LuisBasisLL2.csv');
AllBras = Basis;
AllKets = AllBras';

NMat = csvread('NMat.csv');
absLMat = csvread('absLMat.csv');
littlenMat = csvread('littlenMat.csv');

UMat = csvread('MyUMat.csv');
LMat = csvread('LMat.csv');
VMatnoA = csvread('VMatnoA.csv');

ss=0;
asd=0;
 for Omega = 0.8285  
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
% figure
% plot(OmegaIncreaseMatrix,OmegaMatRestricted);
%  xlabel('Omega')
%  ylabel('<E>')
 end
Dordered = diag(D);

[HamMinVal,HamMinPos] = min(Dordered(:));

[x,y] = size(Basis);
SPDMmat = zeros(y);
creaVect = zeros(y,1);
annVect = zeros(y,1);
 OrigCoeff  = V(:,HamMinPos); 
%  for Kets = 1:x
%      for k = 1:y
%          for l = 1:y
%              C1 =OrigCoeff(Kets,1);             
%              InitialKet = AllKets(:,Kets);
%              creaVect = creaVect.*0;
%              annVect = annVect.*0;
%              creaVect(k,1) = 1;
%              annVect(l,1) = -1;
%              annConst = sqrt(InitialKet(l,1));
%              annTrans = InitialKet + annVect;
%              creaConst = sqrt(annTrans(k,1)+1);
%              FinalKet = annTrans + creaVect;           
%              for Bras = 1:x
%                  Bra = AllBras(Bras,:);
%                  C2 =OrigCoeff(Bras,1);
%                  if Bra == FinalKet'
%                      SPDMmat(k,l) = SPDMmat(k,l) + (creaConst*annConst*C1*C2);
%                  end
%              end
%          end
%      end
%  end
 
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
             if any(annTrans < 0)       
                 continue
             end
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

Vects = [LargestVect SecLargeVect]

TempLarge(:,1) = LargestVect(:);
TempLarge(:,2) = abs(LargestVect(:));
TempSec(:,1) = SecLargeVect(:);
TempSec(:,2) = abs(SecLargeVect(:));

[LmaxVal,LMaxPos] = max(TempLarge(:,2));
TempLarge(LMaxPos,1) = 0;
TempLarge(LMaxPos,2) = 0;
[LSmaxVal,LSmaxPos] = max(TempLarge(:,2));
TempLarge(LSmaxPos,1) = 0;
TempLarge(LSmaxPos,2) = 0;

[SmaxVal,SMaxPos] = max(TempSec(:,2));
TempSec(SMaxPos,1) = 0;
TempSec(SMaxPos,2) = 0;
[SSmaxVal,SSmaxPos] = max(TempSec(:,2));
TempSec(SSmaxPos,1) = 0;
TempSec(SSmaxPos,2) = 0;






toc