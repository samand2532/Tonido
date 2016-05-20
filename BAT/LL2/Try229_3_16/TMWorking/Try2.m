%%% Fresh Start of BatExpansion, using trinomial for GS, Binomial for Ex1

clc; clear all;
tic

KPos = [ 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22; %Matlab position
         0 0 0 0 0 0 0 0 0 0  0  1  1  1  1  1  1  1  1  1  1  1 ; %n value
        -1 0 1 2 3 4 5 6 7 8  9 -1  0  1  2  3  4  5  6  7  8  9]; %Mt value
format long
A = 0.03;
g = 0.4;

Nt = 6;

OmegaLoop = 0;
%for Omega = 0.9:0.000001:0.96
%Omega = 0.9375% GOOD FOR g = 0.4
%Omega = 0.828 g=1 bat????
%Omega = 0.936
Omega = 0.933

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
Ex1VectTemp(Ex1Pos,2) = 0;
[Ex2Val,Ex2Pos] = max(Ex1VectTemp(:,2));
disp('Ex1 Val and Pos')
Ex1Data = [Ex1Val  Ex1Pos; Ex2Val  Ex2Pos]
Ex1Fid = Ex1Val.^2 + Ex2Val.^2

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

[x,y] = size(TriMat);
aa = 0;
for BiCount = 1:x
    N1 = TriMat(BiCount,1);
    N2 = TriMat(BiCount,2);
    i = TriMat(BiCount,3);
    j = TriMat(BiCount,4);
    k = TriMat(BiCount,5);
    
    for p = 0:N2
        aa = aa + 1;
        q = N2-p;
        TriBiMat(aa,1) = N1;
        TriBiMat(aa,2) = N2;
        TriBiMat(aa,3) = i;
        TriBiMat(aa,4) = j;
        TriBiMat(aa,5) = k;
        TriBiMat(aa,6) = p;
        TriBiMat(aa,7) = q;
    end
end

%%%% This section populates the ket in the correct positions with the
%%%% correct values
AllPossKetsN1N2 = zeros(x,24);
ConstAllPosKets = zeros(x,23);
[x,y] = size(TriBiMat);

for KetPop = 1:x;
    N1 = TriBiMat(KetPop,1);
    N2 = TriBiMat(KetPop,2);
    iVal = TriBiMat(KetPop,3);
    jVal = TriBiMat(KetPop,4);
    kVal = TriBiMat(KetPop,5);
    pVal = TriBiMat(KetPop,6);
    qVal = TriBiMat(KetPop,7);
    
    AllPossKetsN1N2(KetPop,1) = N1;
    AllPossKetsN1N2(KetPop,2) = N2;
    AllPossKetsN1N2(KetPop,Gs1Pos+2) = iVal;
    AllPossKetsN1N2(KetPop,Gs2Pos+2) = jVal;
    AllPossKetsN1N2(KetPop,Gs3Pos+2) = kVal;
    AllPossKetsN1N2(KetPop,Ex1Pos+2) = pVal;
    AllPossKetsN1N2(KetPop,Ex2Pos+2) = qVal;
    
    fN1 = factorial(N1);
    fN2 = factorial(N2);
    fi = factorial(iVal);
    fj = factorial(jVal);
    fk = factorial(kVal);
    fq = factorial(qVal);
    fp = factorial(pVal);
    
    ConstAllPosKets(KetPop,1) = (1/(sqrt(fN1*fN2))) * (fN1 / (fi*fj*fk))...
        * (Gs1Val^iVal)*(Gs2Val^jVal)*(Gs3Val^kVal)*sqrt(fi)*sqrt(fj)*sqrt(fk)...
        * nchoosek(N2,pVal)*(Ex1Val^pVal)*(Ex2Val^qVal)*sqrt(fp)*sqrt(fq);
    
    %%% JUST A CHECK TO MAKE SURE CONST IS CORRECT, BOTH AGREE
    %     Val(KetPop,1) = (1/(sqrt(fN1*fN2))) * (fN1 / (fi*fj*fk))...
    %          * (Gs1Val^iVal)*(Gs2Val^jVal)*(Gs3Val^kVal)*sqrt(fi)*sqrt(fj)*sqrt(fk)...
    %          * nchoosek(N2,pVal)*(Ex1Val^pVal)*(Ex2Val^qVal)*sqrt(fp)*sqrt(fq);
    %
    %     ConstAllPosKets(KetPop,1) = (1/(sqrt(fN1*fN2)))*(fN1/(fi*fj*fk))...
    %         *(Gs1Val^iVal)*(Gs2Val^jVal)*(Gs3Val^kVal)*sqrt(fi)*sqrt(fj)*sqrt(fk)...
    %         *nchoosek(N2,pVal)*(Ex1Val^pVal)*(Ex2Val^qVal)*sqrt(fp)*sqrt(fq);
    %     Val(KetPop,2) = (1/(sqrt(fN1*fN2)))*(fN1/(fi*fj*fk))...
    %         *(Gs1Val^iVal)*(Gs2Val^jVal)*(Gs3Val^kVal)*sqrt(fi)*sqrt(fj)*sqrt(fk)...
    %         *nchoosek(N2,pVal)*(Ex1Val^pVal)*(Ex2Val^qVal)*sqrt(fp)*sqrt(fq);
    
    
    
    ConstAllPosKets(KetPop,2:end) = AllPossKetsN1N2(KetPop,3:end);
    
end
%%%% Refines ConstAllPosKets and AllPossKetsN1N2 to exclude daft states, ie
%%%% 6 * -1
aa=0;
for RefineCount = 1:x
    if dot(AllPossKetsN1N2(RefineCount,3:end),KPos(3,:)) < 9
        aa = aa + 1;
        RAllPossKetsN1N2(aa,:) = AllPossKetsN1N2(RefineCount,:);
        RConstAllPosKets(aa,:) = ConstAllPosKets(RefineCount,:);
        %%% RMATNAME for REFINED VERSION OF MATNAME
    end
end
[x,y] = size(RAllPossKetsN1N2);

Ket60Mat=0;
Ket42Mat=0;
Ket24Mat=0;
Ket06Mat=0;
aa=0;bb=0;cc=0;dd=0;
for TMSortCount = 1:x
    N1 = RAllPossKetsN1N2(TMSortCount,1);
    N2 = RAllPossKetsN1N2(TMSortCount,2);
    
    if (N1 == 6) && (N2 == 0)
        for BasisCount = 1:length(Basis)
            if Basis(BasisCount,:) == RConstAllPosKets(TMSortCount,2:end)
                
                aa=aa+1;
                Ket60Mat(aa,1) = RConstAllPosKets(TMSortCount,1)*OrigCoeff(BasisCount,1);
            end
        end
    end
    
    if (N1 == 4) && (N2 == 2)
        for BasisCount = 1:length(Basis)
            if Basis(BasisCount,:) == RConstAllPosKets(TMSortCount,2:end)
                bb=bb+1;
                Ket42Mat(bb,1) = RConstAllPosKets(TMSortCount,1)*OrigCoeff(BasisCount,1);
            end
        end
        
    end
    
    if (N1 == 2) && (N2 == 4)
        for BasisCount = 1:length(Basis)
            if Basis(BasisCount,:) == RConstAllPosKets(TMSortCount,2:end)
                cc=cc+1;
                Ket24Mat(cc,1) = RConstAllPosKets(TMSortCount,1)*OrigCoeff(BasisCount,1);
            end
        end
    end
    
    if (N1 == 0) && (N2 == 6)
        for BasisCount = 1:length(Basis)
            if Basis(BasisCount,:) == RConstAllPosKets(TMSortCount,2:end)
                dd=dd+1;
                Ket06Mat(dd,1) = RConstAllPosKets(TMSortCount,1)*OrigCoeff(BasisCount,1);
            end
        end
    end
    
    if (N1 == 1) && (N2 == 5)
        for BasisCount = 1:length(Basis)
            if Basis(BasisCount,:) == RConstAllPosKets(TMSortCount,2:end)
                dd=dd+1;
                Ket15Mat(dd,1) = RConstAllPosKets(TMSortCount,1)*OrigCoeff(BasisCount,1);
                disp('1,5');
            end
        end
    end
    
    if (N1 == 3) && (N2 == 3)
        for BasisCount = 1:length(Basis)
            if Basis(BasisCount,:) == RConstAllPosKets(TMSortCount,2:end)
                dd=dd+1;
                Ket3Mat(dd,1) = RConstAllPosKets(TMSortCount,1)*OrigCoeff(BasisCount,1);
                disp('3,3');
            end
        end
    end
    
     if (N1 == 5) && (N2 == 1)
        for BasisCount = 1:length(Basis)
            if Basis(BasisCount,:) == RConstAllPosKets(TMSortCount,2:end)
                dd=dd+1;
                Ket3Mat(dd,1) = RConstAllPosKets(TMSortCount,1)*OrigCoeff(BasisCount,1);
                disp('5,1');
            end
        end
    end
    
    
end
Ket60Mat;
Ket42Mat;
Ket24Mat;
Ket06Mat;

wave = [OrigCoeff Basis];
K60Tot = sum(Ket60Mat.^2)
K42Tot = sum(Ket42Mat.^2)
K24Tot = sum(Ket24Mat.^2)
K06Tot = sum(Ket06Mat.^2)

bar([K60Tot K42Tot K24Tot K06Tot])
TotFid = (K60Tot + K42Tot + K24Tot + K06Tot)

OmegaLoop = OmegaLoop + 1;

TempMat(OmegaLoop,1) = Omega;
TempMat(OmegaLoop,2) = K60Tot;
TempMat(OmegaLoop,3) = K42Tot;
TempMat(OmegaLoop,4) = K24Tot;
TempMat(OmegaLoop,5) = K06Tot;
TempMat(OmegaLoop,6) = K60Tot - K06Tot;



toc