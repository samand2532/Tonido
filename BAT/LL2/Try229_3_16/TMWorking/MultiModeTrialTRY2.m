%%%% Multimode trial part 2


clc; clear all;
tic

KPos = [ 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22; %Matlab position
    0 0 0 0 0 0 0 0 0 0  0  1  1  1  1  1  1  1  1  1  1  1 ; %n value
    -1 0 1 2 3 4 5 6 7 8  9 -1  0  1  2  3  4  5  6  7  8  9]; %Mt value
format long;
A = 0.015;
g = 0.6;

Nt = 6;
%%%% Omega = 0.9375722836016 %best so far for g = 0.4
%%%% Omega = somewhere between 0.82856708053814 for g = 1
%%%% Omega = 0.91384786908 for g = 0.5
%%%% Omega = 0.8936197 for g = 0.6
%%%% Omega = 0.86697961104 for g = 0.75
%%%% g = 0.2 Cannot find  an omega that causes Bat state, Omega = 0.99999
%%%% isnt high enough.
%%%% Omega = 0.858808166          g = 0.8   
%%%% Omega = 0.84326424262349    g = 0.9

%%%%%%%%% g Vals Calc g = (0.2??),0.4,0.5,0.6,0.75,0.8,1
start = 0.1
step = 0.1
last = 1.5

 for Count = 1:10
     Count
for Omega = start:step:last
Basis = csvread('EvenBasisLL2.csv');
AllBras = Basis;
AllKets = AllBras';
NMat = csvread('NMat.csv');
absLMat = csvread('absLMat.csv');
littlenMat = csvread('littlenMat.csv');
UMat = csvread('MyUMat.csv');
LMat = csvread('LMat.csv');
VMatnoA = csvread('VMatnoA.csv');

NAll = NMat;
LAll = LMat;
absL = absLMat;

UAll = g.*UMat;
VAll = VMatnoA.*A;

littlenAll = 2.*littlenMat;
Total = littlenAll+(absLMat-(Omega*LAll))+(NAll)+UAll+VAll;
Eig = eig(round(Total,15));

[V,D] = eig(Total);

Dordered = diag(D);

[HamMinVal,HamMinPos] = min(Dordered(:));

[x,y] = size(Basis);
SPDMmat = zeros(y);
creaVect = zeros(y,1);
annVect = zeros(y,1);

OrigCoeff  = V(:,HamMinPos);

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
DiageigVal = diag(eigValSPDM);
[maxValVal,maxValPos] = max(DiageigVal);
VectLarge = eigVectSPDM(:,maxValPos);
DiageigVal(maxValPos,1) = 0;
[SecmaxVal,secMaxPos] = max(DiageigVal);
SecVect = eigVectSPDM(:,secMaxPos);
Vects = [ VectLarge SecVect]

if abs(VectLarge(3,1)) > abs(VectLarge(2,1))
    disp('PAST')
    Omega
    start = Omega - step;
    step = step / 10;
    step
    break
end
end
 end
toc