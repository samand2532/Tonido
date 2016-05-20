%%%%% New expansion

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

for Omega = 0.8258
    %0.818567080538138;
   
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
%% EndSPDM
diageigValSPDM = diag(eigValSPDM);
[maxeigValSPDM,maxPosSPDM] = max(diageigValSPDM);
diageigValSPDM(maxPosSPDM,1) = 0;
[secmaxeigVectSPDM,secmaxPosSPDM] = max(diageigValSPDM);

[maxeigValSPDM,maxPosSPDM];
[secmaxeigVectSPDM,secmaxPosSPDM];

LargestVect = eigVectSPDM(:,maxPosSPDM);
SecLargeVect = eigVectSPDM(:,secmaxPosSPDM);
Vects = [LargestVect SecLargeVect]
LargestVectTemp(:,1) = LargestVect(:);
LargestVectTemp(:,2) = abs(LargestVect(:));
SecLargeVectTemp(:,1) = SecLargeVect(:);
SecLargeVectTemp(:,2) = abs(SecLargeVect(:));

LargestVectTemp;
[maxV1LVal,maxV1Lrow] = max(LargestVectTemp(:,2));
maxV1LVal = LargestVectTemp(maxV1Lrow,1);
LargestVectTemp2 = LargestVectTemp;
LargestVectTemp2(maxV1Lrow,1) = 0;
LargestVectTemp2(maxV1Lrow,2) = 0;
LargestVectTemp2;
[maxV1SVal,maxV1Srow] = max(LargestVectTemp2(:,2));
maxV1SVal = LargestVectTemp2(maxV1Srow,1);

SecLargeVectTemp;
[maxV2LVal,maxV2Lrow] = max(SecLargeVectTemp(:,2));
maxV2LVal = SecLargeVectTemp(maxV2Lrow,1);
SecLargeVectTemp2 = SecLargeVectTemp;
SecLargeVectTemp2(maxV2Lrow,1) = 0;
SecLargeVectTemp2(maxV2Lrow,2) = 0;
SecLargeVectTemp2;
[maxV2SVal,maxV2Srow] = max(SecLargeVectTemp2(:,2));
maxV2SVal = SecLargeVectTemp(maxV2Srow,1);

CL1 = maxV1LVal; CL2 = maxV1SVal; CS1 = maxV2LVal; CS2 = maxV2SVal;
CL1Pos = maxV1Lrow; CL2Pos = maxV1Srow; CS1Pos = maxV2Lrow; CS2Pos = maxV2Srow;
aa=0;
Tab = zeros(100,23);
for N1 = 0:6
    N2 = 6 - N1;
        N1 ;
        N2;
        
    Norm = 1/sqrt((factorial(N1))*factorial(N2));
    
    for p = 0:N1
        for q = 0:N2          
            p;
            q;
            const = (1/((sqrt(factorial(N1))*sqrt(factorial(N2))))) * nchoosek(N1,p) * nchoosek(N2,q) * (CL1^q)*(CL2^p)*...
                (CS1^(N2-q))*(CS2^(N1-p))*sqrt(factorial(p))*sqrt(factorial(q))*...
                sqrt(factorial(N1-p))*sqrt(factorial(N2-q));            
            aa=aa+1;
            Tab(aa,1) = const;
            Tab(aa,CL1Pos+1) = q;
            Tab(aa,CS1Pos+1) = p;
            Tab(aa,CL2Pos+1) = (N2-q);
            Tab(aa,CS2Pos+1) = (N1-p);       
        end
    end        
end

%% This section generates the basis that corresponds with each TM bar position,
N60bar = zeros(1,22);
N42bar = zeros(1,22);
N24bar = zeros(1,22);
N06bar = zeros(1,22);
aa = 0;
for N1 = 6
    N2 = 6 - N1;
        N1 ;
        N2;   
    for p = 0:N1
        for q = 0:N2                            
            aa=aa+1;            
            N60bar(aa,CL1Pos+1) = q;
            N60bar(aa,CS1Pos+1) = p;
            N60bar(aa,CL2Pos+1) = (N2-q);
            N60bar(aa,CS2Pos+1) = (N1-p);       
        end
    end        
end

aa = 0;
for N1 = 4
    N2 = 6 - N1;
        N1 ;
        N2;   
    for p = 0:N1
        for q = 0:N2                            
            aa=aa+1;            
            N42bar(aa,CL1Pos+1) = q;
            N42bar(aa,CS1Pos+1) = p;
            N42bar(aa,CL2Pos+1) = (N2-q);
            N42bar(aa,CS2Pos+1) = (N1-p);       
        end
    end        
end

aa = 0;
for N1 = 2
    N2 = 6 - N1;
        N1 ;
        N2;   
    for p = 0:N1
        for q = 0:N2                            
            aa=aa+1;            
            N24bar(aa,CL1Pos+1) = q;
            N24bar(aa,CS1Pos+1) = p;
            N24bar(aa,CL2Pos+1) = (N2-q);
            N24bar(aa,CS2Pos+1) = (N1-p);       
        end
    end        
end 

aa = 0;
for N1 = 0
    N2 = 6 - N1;
        N1 ;
        N2;   
    for p = 0:N1
        for q = 0:N2                            
            aa=aa+1;            
            N06bar(aa,CL1Pos+1) = q;
            N06bar(aa,CS1Pos+1) = p;
            N06bar(aa,CL2Pos+1) = (N2-q);
            N06bar(aa,CS2Pos+1) = (N1-p);       
        end
    end        
end

%% This section loops Tab (created values using C1 and consts etc) over each of the
%% bar basis to find overlaps
[xTab,yTab] = size(Tab);
[x60,y60] = size(N60bar);
[x42,y42] = size(N42bar);
[x24,y24] = size(N24bar);
[x06,y06] = size(N06bar);
aa=0;bb=0;cc=0;dd=0;ee=0;
for CountA = 1:xTab
    for CountB = 1:x60
        if Tab(CountA,2:6) == N60bar(CountB,2:6)
            for CountW = 1:length(Wavefunction)
                if N60bar(CountB,2:6) == Wavefunction(CountW,2:6)
                    aa=aa+1;
                    N60barVal(aa,1) = Wavefunction(CountW,1) * Tab(CountA,1);
                end
            end
        end
    end
    for CountC = 1:x42
        if Tab(CountA,2:6) == N42bar(CountC,2:6)
            for CountW = 1:length(Wavefunction)
                if N42bar(CountC,2:6) == Wavefunction(CountW,2:6)
                    bb=bb+1;
                    N42barVal(bb,1) = Wavefunction(CountW,1) * Tab(CountA,1);
                end
            end
        end
    end
    for CountD = 1:x24
        if Tab(CountA,2:6) == N24bar(CountD,2:6)
            for CountW = 1:length(Wavefunction)
                if N24bar(CountD,2:6) == Wavefunction(CountW,2:6)
                    cc=cc+1;
                    N24barVal(cc,1) = Wavefunction(CountW,1) * Tab(CountA,1);
                end
            end
        end
    end
    for CountE = 1:x06
        if Tab(CountA,2:6) == N06bar(CountE,2:6)
            for CountW = 1:length(Wavefunction)
                if N06bar(CountE,2:6) == Wavefunction(CountW,2:6)
                    dd=dd+1;
                    N06barVal(dd,1) = Wavefunction(CountW,1) * Tab(CountA,1);
                end
            end
        end
    end
end


bar([sum(N60barVal.^2) sum(N42barVal.^2) sum(N24barVal.^2) sum(N06barVal.^2) ])








toc
































































