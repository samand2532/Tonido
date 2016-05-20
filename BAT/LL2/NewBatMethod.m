%%%% NewBat calculation, using selection to create bat states - so instead
%%%% of |6,0> after OmegaC showing as m=1, it collects these and forces
%%%% them into the |0,6> column.


%%% SPDMmk2 and Expansion Combined

clc; clear all;
tic

KPos = [ 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22; %Matlab position
    0 0 0 0 0 0 0 0 0 0  0   1  1  1  1  1  1  1  1  1  1  1 ; %n value
    -1 0 1 2 3 4 5 6 7 8  9  -1  0  1  2  3  4  5  6  7  8  9]; %Mt value

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
for Omega = 0.9;%0.001:0.9
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

LFidelity = C0Valmax.^2 + C2Valmax.^2 + C4Valmax.^2
SFidelity = C1Valmax.^2 + C3Valmax.^2

%%% L1 = largest of Gs.
%%% L2 = sec large of gs
%%% L3 = third...
%%% S1 = largest small
L1 = C0Valmax;
L1Pos = C0Pos;
L2 = C2Valmax;
L2Pos = C2Pos;
L3 = C4Valmax;
L3Pos = C4Pos;
S1 = C1Valmax;
S1Pos = C1Pos;

%%% expnsion of tri part (L1,L2,L3)
aa = 0;
for i = 0:6
    for j = 0:6
        for k = 0:6
            
            if ((i+j+k) < 7) && ( k < 2)
                aa=aa+1;
                Matijk(aa,1) = (i+j+k);
                Matijk(aa,2) = 6 - (i+j+k);
                Matijk(aa,3) = i;
                Matijk(aa,4) = j;
                Matijk(aa,5) = k;
            end
        end
    end
end

AllBasis = zeros(1,22);
AllBasisTemp = zeros(1,22);
[x,y] = size(Matijk);
AllBasisCount = 0;
N60S=0;
N51S=0;
N42S=0;
N33S=0;
N24S=0;
N15S=0;
N06S=0;
aa=0;
for CountBasis = 1:x
    AllBasisCount = AllBasisCount + 1;
    i = Matijk(CountBasis,3);
    j = Matijk(CountBasis,4);
    k = Matijk(CountBasis,5);
    N1 = i+j+k;
    N2 = N-(i+j+k);
    
    AllBasisTemp(AllBasisCount,L1Pos) =i;
    AllBasisTemp(AllBasisCount,L2Pos) =j;
    AllBasisTemp(AllBasisCount,L3Pos) =k;
    AllBasisTemp(AllBasisCount,S1Pos) = N2;
    
    if dot(AllBasisTemp(CountBasis,:),KPos(3,:)) < 9
        aa = aa + 1;
        AllBasis(aa,:) = AllBasisTemp(CountBasis,:);
    end
end
[x,y] = size(AllBasis);
CountBasis = 0;
for CountA = 1:x
    CountBasis = CountBasis + 1;
    N2 = AllBasis(CountBasis,S1Pos);
    if N2 == 0
        N60S = N60S + 1;
        N60SpecBasis(N60S,:) = AllBasis(CountBasis,:);
    elseif N2 == 1
        N51S = N51S+1;
        N51SpecBasis(N51S,:) = AllBasis(CountBasis,:);
    elseif N2 == 2
        N42S = N42S+1;
        N42SpecBasis(N42S,:) = AllBasis(CountBasis,:);
    elseif N2 == 3
        N33S = N33S+1;
        N33SpecBasis(N33S,:) = AllBasis(CountBasis,:);
    elseif N2 == 4
        N24S = N24S+1;
        N24SpecBasis(N24S,:) = AllBasis(CountBasis,:);
    elseif N2 == 5
        N15S = N15S+1;
        N15SpecBasis(N15S,:) = AllBasis(CountBasis,:);
    elseif N2 == 6
        N06S = N06S+1;
        N06SpecBasis(N06S,:) = AllBasis(CountBasis,:);
    end
    
    factN1 = factorial(N1);
    factN2 = factorial(N2);
    facti = factorial(i);
    factj = factorial(j);
    factk = factorial(k);
    
    ConstVal(CountBasis,1) = (1/((sqrt(factN1))*(sqrt(factN2))))...
        *((factN1)/(facti*factj*factk))*(L1^i)*(L2^j)*(L3^k)...
        *sqrt(facti)*sqrt(factj)*sqrt(factk)*(sqrt(factorial(N2)))*(S1^N2)
end

%%% N60 Ket
[x,y] = size(AllBasis);
[x60,y60] = size(N60SpecBasis);
[x51,y51] = size(N51SpecBasis);
[x42,y42] = size(N42SpecBasis);
[x33,y33] = size(N33SpecBasis);
[x24,y24] = size(N24SpecBasis);
[x15,y15] = size(N15SpecBasis);
[x06,y06] = size(N06SpecBasis);

aa=0;bb=0;cc=0;dd=0;ee=0;ff=0;
for AllCount = 1:x
    for N60Count = 1:x60
        if (N60SpecBasis(N60Count,:) == AllBasis(AllCount,:))
            if  (dot((N60SpecBasis(N60Count,:)),KPos(3,:)) < 9)
                aa =aa+1;
                Ket60Vals(aa,1) = ConstVal(AllCount,1);
                N60SpecBasis(N60Count,:)
                AllBasis(AllCount,:)
            end
        end
    end
    for N51Count = 1:x51
        if N51SpecBasis(N51Count,:) == AllBasis(AllCount,:)
            if (dot((N51SpecBasis(N51Count,:)),KPos(3,:)) < 9)
                bb = bb+1;
                Ket51Vals(bb,1) = ConstVal(AllCount,1);
            end
        end
    end
    for N42Count = 1:x42
        if N42SpecBasis(N42Count,:) == AllBasis(AllCount,:)
            if (dot((N42SpecBasis(N42Count,:)),KPos(3,:)) < 9)
                cc = cc+1;
                Ket42Vals(cc,1) = ConstVal(AllCount,1);
            end
        end
    end
    for N33Count = 1:x33
        if N33SpecBasis(N33Count,:) == AllBasis(AllCount,:)
            if (dot((N33SpecBasis(N33Count,:)),KPos(3,:)) < 9)
                dd = dd + 1;
                Ket33Vals(dd,1) = ConstVal(AllCount,1);
            end
        end
    end
    for N24Count = 1:x24
        if N24SpecBasis(N24Count,:) == AllBasis(AllCount,:)
            if  (dot((N24SpecBasis(N24Count,:)),KPos(3,:)) < 9)
                dd = dd+1;
                Ket24Vals(dd,1) = ConstVal(AllCount,1);
            end
        end
    end
    for N15Count = 1:x15
        if N15SpecBasis(N15Count,:) == AllBasis(AllCount,:)
            if (dot((N15SpecBasis(N15Count,:)),KPos(3,:)) < 9)
                ee = ee+1;
                Ket15Vals(ee,1) = ConstVal(AllCount,1);
            end
        end
    end
    for N06Count = 1:x06
        if N06SpecBasis(N06Count,:) == AllBasis(AllCount,:)
            if (dot((N06SpecBasis(N06Count,:)),KPos(3,:)) < 9)
                ff = ff+1;
                Ket06Vals(ff,1) = ConstVal(AllCount,1);
            end
        end
    end
end
KetN60Val = sum(Ket60Vals(:).^2)
KetN51Val = sum(Ket51Vals(:).^2)
KetN42Val = sum(Ket42Vals(:).^2)
KetN33Val = sum(Ket33Vals(:).^2)
KetN24Val = sum(Ket24Vals(:).^2)
KetN15Val = sum(Ket15Vals(:).^2)
KetN06Val = sum(Ket06Vals(:).^2)

FidelityOverlap = KetN60Val + KetN51Val + KetN42Val + KetN33Val...
    + KetN24Val + KetN15Val + KetN06Val

toc






