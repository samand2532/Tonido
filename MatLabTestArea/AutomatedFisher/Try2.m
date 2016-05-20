tic
clear all; clc;
DeleteMatrix = 1; % 1 = yes
N = 6;Nt = 6; g =1; LLVal = 1;% For the deltainteractionfunction, LLL = 1, LL2incLLL = 2
A = 0.03; 
MomValue = [ -1 0 0 1 1 2 2 3 3 4 4 5 5 6 6 7 7 8 8]'; 
if LLVal == 1
AllBras = csvread('StatesLLL.csv');
elseif LLVavl ==2
    AllBras = csvread('StatesLLL_and_LL2.csv');
end
AllKets = AllBras';

LengthAllBras = length(AllBras);
LengthVect = length(AllBras(1,:));
NMatLLL(N, LengthAllBras); %Function
LMatLLL(LLVal); % Function
DeltaInteractionVariableLL(LLVal, LengthVect);%Function
InteractionTermLLL(g, LengthVect, LLVal);%Function
PotentialLLL(LengthVect, LLVal); %Function
ss=0;
bb = 0;
for Omega = 0.782
    ss=ss+1   
    NAll = csvread('NMatLLLFinal.csv');
    UAll = csvread('UMatFinalLLL.csv');
    LAll = (csvread('LMatLLL.csv')).*(1-Omega);
    VNoA = csvread('VMatNoA_LLL.csv');
    VAll = VNoA.*0.5.*A;
    
    Total = NAll+UAll+LAll+VAll;
    Eig = eig(Total);
    [V, D] = eig(Total);
    EigStateMat(1,ss) = Omega;
    EigStateMat(2:40,ss) = V(:,1);  
    SPDMMat = zeros(3,4010);
FinalSortMat = zeros(2,1);

bb = bb + 1;
AllBras = csvread('StatesLLL.csv');
AllKets = AllBras';
V = EigStateMat(2:40,bb);
Omega = EigStateMat(1,bb);
eigenstates = V(:,1);
format long
eigMat = zeros(9);

creaVect = zeros(19,1);
annVect = zeros(19,1);

for Kets = 1:39;
for k = [2 4 6 8 10 12 14 16 18];
for l = [2 4 6 8 10 12 14 16 18];
    C1 = eigenstates(Kets,1);
    InitialKet = AllKets(:,Kets);
    creaVect = creaVect.*0;
    annVect = annVect.*0;
    creaVect(k,1) = 1;
    annVect(l,1) = -1;
    
    annConst = sqrt(InitialKet(l,1));
    annTrans = InitialKet + annVect;
    creaConst = sqrt(annTrans(k,1)+1);
    FinalKet = annTrans + creaVect;
a=0;
    for a = 1:39
        
        Bra = AllBras(a,:);
        C2 = eigenstates(a,1);
        if Bra' == FinalKet
            eigMat(k/2,l/2) = eigMat(k/2,l/2) + (creaConst*annConst*C1*C2);
        end
    end
end
end
end

[V1 D] = eig(eigMat);

[M I] = max(D(:));
[row col] = ind2sub(size(D),I);

ATemp = V1(:,col);
ATempAbs = abs(ATemp);

[MR IR] = max(ATempAbs); %%% FInds max and position of Max
[Rrow Rcol] = ind2sub(size(ATempAbs),IR);
BinLarge = ATemp(Rrow,Rcol);
ATempAbs(Rrow,Rcol) = 0; %%%% Line wipes out max value to allow 2nd highest to be found

[MS IS] = max(ATempAbs); %%% FInds max and position of Max
[Srow Scol] = ind2sub(size(ATempAbs),IS);
BinSmaller = ATemp(Srow,Scol);

SPDMMat(1,bb) = Omega;
SPDMMat(2,bb) = BinLarge;
SPDMMat(3,bb) = BinSmaller;

%%%%%%%%%% WORKING UPTO HERE

x = BinLarge;
x1 = BinSmaller;
Eig = EigStateMat(2:40,bb);



    
    k=0;
    KetN0=zeros(1,10);
for N0 = 0
    for p = 0:N0
        k=k+1;
        
        Nk = nchoosek(N0,p);
        Norm = 1/((sqrt(factorial(N0))).*(sqrt(factorial(Nt-N0))));
        
        
        KetN0(k,1) = Nk.*Norm.*(x^p).*(x1^(N0-p)).*(sqrt(factorial(p)))...
            .*sqrt(factorial(N0-p)).*sqrt(factorial(Nt-N0));
        KetN0(k,2) = p;
        KetN0(k,3) = 6-N0;
        KetN0(k,4) = N0-p;
    end
end

k=0;
KetN1 = zeros(2,10);
for N1 = 1
    for p = 0:N1
        k=k+1;

         Nk = nchoosek(N1,p)
         SumTerm = Nk.*((x)^(N1-p)).*((x1).^p)
         Norm = 1/((sqrt(factorial(N1))).*(sqrt(factorial(6-N1))))
         TotConst = SumTerm.*Norm.*Nk


        KetN1(k,1) = Nk.*SumTerm.*Norm*(sqrt(factorial(p))*sqrt(factorial(6-N1))*sqrt(factorial(N1-p)));
        KetN1(k,2) = p;
        KetN1(k,3) = 6-N1;
        KetN1(k,4) = N1-p;
    end
end

k=0;
KetN2 = zeros(3,10);
for N2 = 2
    for p = 0:N2
        k=k+1;
        
        Nk = nchoosek(N2,p);
        Norm = 1/((sqrt(factorial(N2))).*(sqrt(factorial(Nt-N2))));
        
        
        KetN2(k,1) = Nk.*Norm.*(x^p).*(x1^(N2-p)).*(sqrt(factorial(p)))...
            .*sqrt(factorial(N2-p)).*sqrt(factorial(Nt-N2));
        KetN2(k,2) = p;
        KetN2(k,3) = 6-N2;
        KetN2(k,4) = N2-p;
    end
end


k=0;
KetN3 = zeros(4,10);
for N3 = 3
    for p = 0:N3
        k=k+1;

         Nk = nchoosek(N3,p)
         SumTerm = Nk.*((x)^(N3-p)).*((x1).^p)
         Norm = 1/((sqrt(factorial(N3))).*(sqrt(factorial(6-N3))))
         TotConst = SumTerm.*Norm.*Nk


        KetN3(k,1) = Nk.*SumTerm.*Norm*(sqrt(factorial(p))*sqrt(factorial(6-N3))*sqrt(factorial(N3-p)));
        KetN3(k,2) = p;
        KetN3(k,3) = 6-N3;
        KetN3(k,4) = N3-p;
    end
end

k=0;
KetN4 = zeros(5,10);
for N4 = 4
    for p = 0:N4
        k=k+1;
        
        Nk = nchoosek(N4,p);
        Norm = 1/((sqrt(factorial(N4))).*(sqrt(factorial(Nt-N4))));
        
        
        KetN4(k,1) = Nk.*Norm.*(x^p).*(x1^(N4-p)).*(sqrt(factorial(p)))...
            .*sqrt(factorial(N4-p)).*sqrt(factorial(Nt-N4));
        KetN4(k,2) = p;
        KetN4(k,3) = 6-N4;
        KetN4(k,4) = N4-p;
    end
end

k=0;
KetN5 = zeros(6,10);
for N5 = 5
    for p = 0:N5
        k=k+1;

         Nk = nchoosek(N5,p)
         SumTerm = Nk.*((x)^(N5-p)).*((x1).^p)
         Norm = 1/((sqrt(factorial(N5))).*(sqrt(factorial(6-N5))))
         TotConst = SumTerm.*Norm.*Nk


        KetN5(k,1) = Nk.*SumTerm.*Norm*(sqrt(factorial(p))*sqrt(factorial(6-N5))*sqrt(factorial(N5-p)));
        KetN5(k,2) = p;
        KetN5(k,3) = 6-N5;
        KetN5(k,4) = N5-p;
    end
end


k=0;
KetN6 = zeros(7,10);
for N6 = 6
    for p = 0:N6
        k=k+1;
        
        Nk = nchoosek(N6,p);
        Norm = 1/((sqrt(factorial(N6))).*(sqrt(factorial(Nt-N6))));
        
        
        KetN6(k,1) = Nk.*Norm.*(x^p).*(x1^(N6-p)).*(sqrt(factorial(p)))...
            .*sqrt(factorial(N6-p)).*sqrt(factorial(Nt-N6));
        KetN6(k,2) = p;
        KetN6(k,3) = 6-N6;
        KetN6(k,4) = N6-p;
    end
end

BasisChangeKet = [KetN0;KetN2;KetN4;KetN6];

%%%%%Above generates the basis transformation from |N_1,N_2> to the
%%%%%momentum basis used initially, next part finds fidelity with riginal
%%%%%wavefunction.


OrigStates = csvread('StatesLLL.csv');
OrigWavefunction = [Eig(:,1) OrigStates];

%%%%%% Ket0 Part
l=0;m=0;
KetN0SaveMat = zeros(10,1);
for CountA = 1:length(OrigWavefunction)
    l=l+1;
    if KetN0(1,2:10) == OrigWavefunction(l,2:10)
        m=m+1;
        KetN0SaveMat(m,1) = KetN0(1,1)*OrigWavefunction(l,1) ;
    end
end

%%%%%%% Ket2 Part
l=0;m=0;
KetN2SaveMat = zeros(10,1);
for CountA = 1:length(OrigWavefunction)
    l=l+1;
    if KetN2(1,2:10) == OrigWavefunction(l,2:10)
        m=m+1;
        KetN2SaveMat(m,1) = KetN2(1,1)*OrigWavefunction(l,1);
    elseif KetN2(2,2:10) == OrigWavefunction(l,2:10)
        m=m+1;
        KetN2SaveMat(m,1) = KetN2(2,1)*OrigWavefunction(l,1);
    elseif KetN2(3,2:10) == OrigWavefunction(l,2:10)
        m=m+1;
        KetN2SaveMat(m,1) = KetN2(3,1)*OrigWavefunction(l,1);
    end
end

%%%%%%% Ket4 Part
l=0;m=0;
KetN4SaveMat = zeros(10,1);
for CountA = 1:length(OrigWavefunction)
    l=l+1;
    if KetN4(1,2:10) == OrigWavefunction(l,2:10)
        m=m+1;
        KetN4SaveMat(m,1) = KetN4(1,1)*OrigWavefunction(l,1);
    elseif KetN4(2,2:10) == OrigWavefunction(l,2:10)
        m=m+1;
        KetN4SaveMat(m,1) = KetN4(2,1)*OrigWavefunction(l,1);
    elseif KetN4(3,2:10) == OrigWavefunction(l,2:10)
        m=m+1;
        KetN4SaveMat(m,1) = KetN4(3,1)*OrigWavefunction(l,1);
    elseif KetN4(4,2:10) == OrigWavefunction(l,2:10)
        m=m+1;
        KetN4SaveMat(m,1) = KetN4(4,1)*OrigWavefunction(l,1);
    elseif KetN4(5,2:10) == OrigWavefunction(l,2:10)
        m=m+1;
        KetN4SaveMat(m,1) = KetN4(5,1)*OrigWavefunction(l,1);
    end
end

%%%%%Ket6 Part
l=0;m=0;
KetN6SaveMat = zeros(10,1);
for CountA = 1:length(OrigWavefunction)
    l=l+1;
    if KetN6(1,2:10) == OrigWavefunction(l,2:10)
        m=m+1;
        KetN6(1,1);
        OrigWavefunction(l,1);
        KetN6SaveMat(m,1) = KetN6(1,1)*OrigWavefunction(l,1);
    elseif KetN6(2,2:10) == OrigWavefunction(l,2:10)
        m=m+1;
        KetN6(2,1);
        OrigWavefunction(l,1);
        KetN6SaveMat(m,1) = KetN6(2,1)*OrigWavefunction(l,1);
    elseif KetN6(3,2:10) == OrigWavefunction(l,2:10)
        m=m+1;
        KetN6(3,1);
        OrigWavefunction(l,1);
        KetN6SaveMat(m,1) = KetN6(3,1)*OrigWavefunction(l,1);
    elseif KetN6(4,2:10) == OrigWavefunction(l,2:10)
        m=m+1;
        KetN6(4,1);
        OrigWavefunction(l,1);
        KetN6SaveMat(m,1) = KetN6(4,1)*OrigWavefunction(l,1);
    elseif KetN6(5,2:10) == OrigWavefunction(l,2:10)
        m=m+1;
        KetN6(5,1);
        OrigWavefunction(l,1);
        KetN6SaveMat(m,1) = KetN6(5,1)*OrigWavefunction(l,1);
    elseif KetN6(6,2:10) == OrigWavefunction(l,2:10)
        m=m+1;
        KetN6(6,1);
        OrigWavefunction(l,1);
        KetN6SaveMat(m,1) = KetN6(6,1)*OrigWavefunction(l,1);
    elseif KetN6(7,2:10) == OrigWavefunction(l,2:10)
        m=m+1;
        KetN6(7,1);
        OrigWavefunction(l,1);
        KetN6SaveMat(m,1) = KetN6(7,1)*OrigWavefunction(l,1);
    end
end

KetN0Overlap = (abs(sum(KetN0SaveMat)))^2;
KetN2Overlap = (abs(sum(KetN2SaveMat)))^2;
KetN4Overlap = (abs(sum(KetN4SaveMat)))^2;
KetN6Overlap = (abs(sum(KetN6SaveMat)))^2;

Overlaps = [KetN0Overlap KetN2Overlap KetN4Overlap KetN6Overlap]
bar(Overlaps)
hold on
    
end


























