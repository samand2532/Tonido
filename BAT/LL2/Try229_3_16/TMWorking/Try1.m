%%% Fresh Start of BatExpansion, using trinomial for GS, only single for
%%% Ex1

clc; clear all;
tic
A = 0.03;
g = 0.4;

Nt = 6;
Omega = 0.933698


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
format short

[eigVectSPDM,eigValSPDM] = eig(SPDMmat);

%%%%% finds 2 largest (+ve) eigenstates
%%% by diaganolising (MATLAB DEFINITION) the eigValSPDM mat, finding the 
%%% corresponding row, then overwriting it with 0 then
%%% finding the next max.

diaEigVal = diag(eigValSPDM);
[maxEigVal,maxEigPos] = max(diaEigVal);
diaEigVal(maxEigPos,1) = 0;
[SecmaxEigVal,SecmaxEigPos] = max(diaEigVal);
GsVect = eigVectSPDM(:,maxEigPos);
Ex1Vect = eigVectSPDM(:,SecmaxEigPos);
GsVectTemp(:,1) = (GsVect(:));
GsVectTemp(:,2) = abs(GsVect(:));
[Gs1Val,Gs1Pos] = max(GsVectTemp(:,2));
Gs1Val = GsVectTemp(Gs1Pos,1);
GsVectTemp(Gs1Pos,:) = 0;
[Gs2Val,Gs2Pos] = max(GsVectTemp(:,2));
Gs2Val = GsVectTemp(Gs2Pos,1);
GsVectTemp(Gs2Pos,:) = 0;
[Gs3Val,Gs3Pos] = max(GsVectTemp(:,2));
Gs3Val = GsVectTemp(Gs3Pos,1);
GsVectTemp(Gs3Pos,:) = 0;
Vects = [GsVect Ex1Vect]
disp('Gs Val and Pos')
GsData = [Gs1Val Gs1Pos;Gs2Val Gs2Pos;Gs3Val Gs3Pos]
GsFid = (Gs1Val.^2)+(Gs2Val.^2)+(Gs3Val.^2)

Ex1VectTemp(:,1) = Ex1Vect(:,1);
Ex1VectTemp(:,2) = abs(Ex1Vect(:,1));
[Ex1Val,Ex1Pos] = max(Ex1VectTemp(:,2));
Ex1Val = Ex1VectTemp(Ex1Pos,1);
disp('Ex1 Val and Pos')
Ex1Data = [Ex1Val  Ex1Pos]
Ex1Fid = Ex1Val.^2

%%%% Trinomail expansion part
aa = 0;
for loopi = 0:6
    for loopj = 0:6
        for loopk = 0:6
            if (loopi+loopj+loopk < 7)
                aa = aa +1;
            TriMat(aa,1) = loopi+loopj+loopk; %N1
            TriMat(aa,2) = 6 - (loopi+loopj+loopk); %N2
            TriMat(aa,3) = loopi;
            TriMat(aa,4) = loopj;
            TriMat(aa,5) = loopk;
            end
        end
    end
end

AllPosKets=zeros(length(TriMat),22);
AllPosKetsN1N2=zeros(length(TriMat),24);
for CountA = 1:length(TriMat)
    N1 = TriMat(CountA,1);
    N2 = TriMat(CountA,2);
    iVal = TriMat(CountA,3);
    jVal = TriMat(CountA,4);
    kVal = TriMat(CountA,5);
    AllPosKets(CountA,Gs1Pos) = iVal;
    AllPosKets(CountA,Gs2Pos) = jVal;
    AllPosKets(CountA,Gs3Pos) = kVal;
    AllPosKets(CountA,Ex1Pos) = N2;
    AllPosKetsN1N2(CountA,1) = N1;
    AllPosKetsN1N2(CountA,2) = N2;
    AllPosKetsN1N2(CountA,Gs1Pos+2) = iVal;
    AllPosKetsN1N2(CountA,Gs2Pos+2) = jVal;
    AllPosKetsN1N2(CountA,Gs3Pos+2) = kVal;
    AllPosKetsN1N2(CountA,Ex1Pos+2) = N2;
end

for CountB = 1:length(TriMat)
    N1 = TriMat(CountB,1);
    N2 = TriMat(CountB,2);
    iVal = TriMat(CountB,3);
    jVal = TriMat(CountB,4);
    kVal = TriMat(CountB,5);
    ConstAllPosKets(CountB,Gs1Pos+1) = iVal;
    ConstAllPosKets(CountB,Gs2Pos+1) = jVal;
    ConstAllPosKets(CountB,Gs3Pos+1) = kVal;
    ConstAllPosKets(CountB,Ex1Pos+1) = N2;
    
    fN1 = factorial(N1);
    fN2 = factorial(N2);
    fi = factorial(iVal);
    fj = factorial(jVal);
    fk = factorial(kVal);
    
    ConstAllPosKets(CountB,1) = (1/((sqrt(fN1))*(sqrt(fN2))))...
        *((fN1)/(fi*fj*fk))*(Gs1Val^iVal)*(Gs2Val^jVal)*(Gs3Val^kVal)...
        *sqrt(fi)*sqrt(fj)*sqrt(fk)...
        *(Ex1Val^N2)*sqrt(fN2);
end

Ket60Mat=0;
Ket42Mat=0;
Ket24Mat=0;
Ket06Mat=0;
aa=0;bb=0;cc=0;dd=0;
for CountC = 1:length(AllPosKetsN1N2)
    N1 = AllPosKetsN1N2(CountC,1);
    N2 = AllPosKetsN1N2(CountC,2);
    
    if (N1 == 6) && (N2 == 0)
        %disp('60');
        
        for BasisCount = 1:length(Basis)
            if Basis(BasisCount,1:22) == AllPosKetsN1N2(CountC,3:24)
                aa=aa+1;
                Ket60Mat(aa,1) = ConstAllPosKets(CountC,1)*OrigCoeff(BasisCount,1);
            end
        end
    end


    if (N1 == 4) && (N2 == 2)
        %disp('42');
        for BasisCount = 1:length(Basis)
            if Basis(BasisCount,1:22) == AllPosKetsN1N2(CountC,3:24)
                bb=bb+1;
                Ket42Mat(bb,1) = ConstAllPosKets(CountC,1)*OrigCoeff(BasisCount,1);
            end
        end
    end
    
    
    
    if (N1 == 2 ) && (N2 == 4)
        %disp('24');
        for BasisCount = 1:length(Basis)
            if Basis(BasisCount,1:22) == AllPosKetsN1N2(CountC,3:24)
                cc=cc+1;
                Ket24Mat(cc,1) = ConstAllPosKets(CountC,1)*OrigCoeff(BasisCount,1);
            end
        end
    end
    
    
    
    if (N1 == 0) && (N2 == 6)
        %disp('06');
        for BasisCount = 1:length(Basis)
            if Basis(BasisCount,1:22) == AllPosKetsN1N2(CountC,3:24)
                dd=dd+1;
                Ket06Mat(dd,1) = ConstAllPosKets(CountC,1)*OrigCoeff(BasisCount,1);
                BasisCount
                CountC
            end
        end
    end





end

P60 = sum(Ket60Mat(:).^2)
P42 = sum(Ket42Mat(:).^2)
P24 = sum(Ket24Mat(:).^2)
P06 = sum(Ket06Mat(:).^2)

%bar([P60 P42 P24 P06])

%TotFid = P60 + P42 + P24 + P06

W = 






























toc