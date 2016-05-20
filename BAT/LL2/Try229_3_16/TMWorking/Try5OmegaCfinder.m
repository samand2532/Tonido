%%%%%% Section selects the number of states to calc depending on if they
%%%%%% exceed a pre-defined tolerance, ie > 0.97



clc; clear all;
tic

KPos = [ 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22; %Matlab position
    0 0 0 0 0 0 0 0 0 0  0  1  1  1  1  1  1  1  1  1  1  1 ; %n value
    -1 0 1 2 3 4 5 6 7 8  9 -1  0  1  2  3  4  5  6  7  8  9]; %Mt value
format long;
A = 0.03;
g = 1;

Nt = 6;
%%%% Omega = 0.9375722836016 %best so far for g = 0.4
%%%% Omega = somewhere between 0.828567080500000 and 0.828567080600000 for g = 1

aaa = 0.828567080500000;
ccc = 0.828567080600000;
bbb = (ccc-aaa)/10;
datestr(now)
for Omega = aaa:bbb:ccc
    Omega
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
%%%%% Section finds 2 largest pops of eigVal, then finds the amkeup of
%%%%% vects, and closes gap towards Omega_c
[maxEVal,maxEPos] = max(diag(eigValSPDM));
TestVect = eigVectSPDM(:,maxEPos)
if abs(TestVect(3,1)) > abs(TestVect(2,1))
    disp('PAST')
    aaa = Omega - bbb
    ccc = Omega
    bbb = (ccc-aaa)/10
    break
end
end
datestr(now)
for Omega = aaa:bbb:ccc
    Omega
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
%%%%% Section finds 2 largest pops of eigVal, then finds the amkeup of
%%%%% vects, and closes gap towards Omega_c
[maxEVal,maxEPos] = max(diag(eigValSPDM));
TestVect = eigVectSPDM(:,maxEPos)
if abs(TestVect(3,1)) > abs(TestVect(2,1))
    disp('PAST')
    aaa = Omega - bbb
    ccc = Omega
    bbb = (ccc-aaa)/10
    break
end
end
datestr(now)
for Omega = aaa:bbb:ccc
    Omega
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
%%%%% Section finds 2 largest pops of eigVal, then finds the amkeup of
%%%%% vects, and closes gap towards Omega_c
[maxEVal,maxEPos] = max(diag(eigValSPDM));
TestVect = eigVectSPDM(:,maxEPos)
if abs(TestVect(3,1)) > abs(TestVect(2,1))
    disp('PAST')
    aaa = Omega - bbb
    ccc = Omega
    bbb = (ccc-aaa)/10
    break
end
end
datestr(now)
for Omega = aaa:bbb:ccc
    Omega
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
%%%%% Section finds 2 largest pops of eigVal, then finds the amkeup of
%%%%% vects, and closes gap towards Omega_c
[maxEVal,maxEPos] = max(diag(eigValSPDM));
TestVect = eigVectSPDM(:,maxEPos)
if abs(TestVect(3,1)) > abs(TestVect(2,1))
    disp('PAST')
    aaa = Omega - bbb
    ccc = Omega
    bbb = (ccc-aaa)/10
    break
end
end
datestr(now)
for Omega = aaa:bbb:ccc
    Omega
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
%%%%% Section finds 2 largest pops of eigVal, then finds the amkeup of
%%%%% vects, and closes gap towards Omega_c
[maxEVal,maxEPos] = max(diag(eigValSPDM));
TestVect = eigVectSPDM(:,maxEPos)
if abs(TestVect(3,1)) > abs(TestVect(2,1))
    disp('PAST')
    aaa = Omega - bbb
    ccc = Omega
    bbb = (ccc-aaa)/10
    break
end
end
datestr(now)
for Omega = aaa:bbb:ccc
    Omega
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
%%%%% Section finds 2 largest pops of eigVal, then finds the amkeup of
%%%%% vects, and closes gap towards Omega_c
[maxEVal,maxEPos] = max(diag(eigValSPDM));
TestVect = eigVectSPDM(:,maxEPos)
if abs(TestVect(3,1)) > abs(TestVect(2,1))
    disp('PAST')
    aaa = Omega - bbb
    ccc = Omega
    bbb = (ccc-aaa)/10
    break
end
end
datestr(now)
for Omega = aaa:bbb:ccc
    Omega
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
%%%%% Section finds 2 largest pops of eigVal, then finds the amkeup of
%%%%% vects, and closes gap towards Omega_c
[maxEVal,maxEPos] = max(diag(eigValSPDM));
TestVect = eigVectSPDM(:,maxEPos)
if abs(TestVect(3,1)) > abs(TestVect(2,1))
    disp('PAST')
    aaa = Omega - bbb
    ccc = Omega
    bbb = (ccc-aaa)/10
    break
end
end
datestr(now)
for Omega = aaa:bbb:ccc
    Omega
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
%%%%% Section finds 2 largest pops of eigVal, then finds the amkeup of
%%%%% vects, and closes gap towards Omega_c
[maxEVal,maxEPos] = max(diag(eigValSPDM));
TestVect = eigVectSPDM(:,maxEPos)
if abs(TestVect(3,1)) > abs(TestVect(2,1))
    disp('PAST')
    aaa = Omega - bbb
    ccc = Omega
    bbb = (ccc-aaa)/10
    break
end
end












toc