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
for Omega = 0.8285670805
    %0.828567080538137;%0.001:0.9
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

%% expansion of tri part / allowed vals of i,j,k
aa = 0;
for i = 0:6
    for j = 0:6
        for k = 0:6
            if ((i+j+k)) <7 && (k < 2)
                
                %%%% i|j|k|N1|N2
            aa = aa + 1;
            Matijk(aa,1) = i;
            Matijk(aa,2) = j;
            Matijk(aa,3) = k;
            Matijk(aa,4) = i+j+k;
            Matijk(aa,5) = N-(i+j+k);
                
            end
        end
    end
end

%%%% expansion of Bi Parrt
bb = 0; cc = 0; dd = 0; ee = 0; ff = 0; gg = 0; hh = 0;
[x,y] = size(Matijk)
for CountA = 1:x
    N1 = Matijk(CountA,4);
    N2 = Matijk(CountA,5);
    i = Matijk(CountA,1);
    j = Matijk(CountA,2);
    k = Matijk(CountA,3);
    factN1 = factorial(N1);
    factN2 = factorial(N2);
    facti = factorial(i);
    factj = factorial(j);
    factk = factorial(k);
    if N2 == 0
        for p = 0:N2
            bb = bb +1;
            Ket60(bb,1) = (1/((sqrt(factN1))*(sqrt(factN2))))...
                *((factN1)/(facti*factj*factk))*nchoosek(N2,p)...
                *(C0Val^i)*(C2Val^j)*(C4Val^k)*(C1Val^p)*(C3Val^(N2-p))...
                *sqrt(facti)*sqrt(factj)*sqrt(factk)*sqrt(factorial(p))*sqrt(factorial(N2-p));
            
            Ket60(bb,C1Pos+1) = p;
            Ket60(bb,C3Pos+1) = (N2-p);
            Ket60(bb,C0Pos+1) = i;
            Ket60(bb,C2Pos+1) = j;
            Ket60(bb,C4Pos+1) = k;
        end
        
        
    end
    
    if N2 == 1
        for p = 0:N2
            cc = cc +1;
            Ket51(cc,1) = (1/((sqrt(factN1))*(sqrt(factN2))))...
                *((factN1)/(facti*factj*factk))*nchoosek(N2,p)...
                *(C0Val^i)*(C2Val^j)*(C4Val^k)*(C1Val^p)*(C3Val^(N2-p))...
                *sqrt(facti)*sqrt(factj)*sqrt(factk)*sqrt(factorial(p))*sqrt(factorial(N2-p));
            
            Ket51(cc,C1Pos+1) = p;
            Ket51(cc,C3Pos+1) = (N2-p);
            Ket51(cc,C0Pos+1) = i;
            Ket51(cc,C2Pos+1) = j;
            Ket51(cc,C4Pos+1) = k;
        end
     
    end
    
    if N2 == 2
        for p = 0:N2
            dd = dd +1;
            Ket42(dd,1) = (1/((sqrt(factN1))*(sqrt(factN2))))...
                *((factN1)/(facti*factj*factk))*nchoosek(N2,p)...
                *(C0Val^i)*(C2Val^j)*(C4Val^k)*(C1Val^p)*(C3Val^(N2-p))...
                *sqrt(facti)*sqrt(factj)*sqrt(factk)*sqrt(factorial(p))*sqrt(factorial(N2-p));
            Ket42(dd,C1Pos+1) = p;
            Ket42(dd,C3Pos+1) = (N2-p);
            Ket42(dd,C0Pos+1) = i;
            Ket42(dd,C2Pos+1) = j;
            Ket42(dd,C4Pos+1) = k;
        end
        
        
    end
    
    if N2 == 3
        for p = 0:N2
            ee = ee +1;
            Ket33(ee,1) = (1/((sqrt(factN1))*(sqrt(factN2))))...
                *((factN1)/(facti*factj*factk))*nchoosek(N2,p)...
                *(C0Val^i)*(C2Val^j)*(C4Val^k)*(C1Val^p)*(C3Val^(N2-p))...
                *sqrt(facti)*sqrt(factj)*sqrt(factk)*sqrt(factorial(p))*sqrt(factorial(N2-p));
            Ket33(ee,C1Pos+1) = p;
            Ket33(ee,C3Pos+1) = (N2-p);
            Ket33(ee,C0Pos+1) = i;
            Ket33(ee,C2Pos+1) = j;
            Ket33(ee,C4Pos+1) = k;
        end
        
        
    end
    
    if N2 == 4
        for p = 0:N2
            ff = ff +1;
            Ket24(ff,1) = (1/((sqrt(factN1))*(sqrt(factN2))))...
                *((factN1)/(facti*factj*factk))*nchoosek(N2,p)...
                *(C0Val^i)*(C2Val^j)*(C4Val^k)*(C1Val^p)*(C3Val^(N2-p))...
                *sqrt(facti)*sqrt(factj)*sqrt(factk)*sqrt(factorial(p))*sqrt(factorial(N2-p));
            Ket24(ff,C1Pos+1) = p;
            Ket24(ff,C3Pos+1) = (N2-p);
            Ket24(ff,C0Pos+1) = i;
            Ket24(ff,C2Pos+1) = j;
            Ket24(ff,C4Pos+1) = k;
        end
        
        
    end
    
    if N2 == 5
        for p = 0:N2
            gg = gg +1;
            Ket15(gg,1) = (1/((sqrt(factN1))*(sqrt(factN2))))...
                *((factN1)/(facti*factj*factk))*nchoosek(N2,p)...
                *(C0Val^i)*(C2Val^j)*(C4Val^k)*(C1Val^p)*(C3Val^(N2-p))...
                *sqrt(facti)*sqrt(factj)*sqrt(factk)*sqrt(factorial(p))*sqrt(factorial(N2-p));
            Ket15(gg,C1Pos+1) = p;
            Ket15(gg,C3Pos+1) = (N2-p);
            Ket15(gg,C0Pos+1) = i;
            Ket15(gg,C2Pos+1) = j;
            Ket15(gg,C4Pos+1) = k;
        end
        
        
    end
    
    if N2 == 6
        for p = 0:N2
            hh = hh +1;
            Ket06(hh,1) = (1/((sqrt(factN1))*(sqrt(factN2))))...
                *((factN1)/(facti*factj*factk))*nchoosek(N2,p)...
                *(C0Val^i)*(C2Val^j)*(C4Val^k)*(C1Val^p)*(C3Val^(N2-p))...
                *sqrt(facti)*sqrt(factj)*sqrt(factk)*sqrt(factorial(p))*sqrt(factorial(N2-p));
            Ket06(hh,C1Pos+1) = p;
            Ket06(hh,C3Pos+1) = (N2-p);
            Ket06(hh,C0Pos+1) = i;
            Ket06(hh,C2Pos+1) = j;
            Ket06(hh,C4Pos+1) = k;
            
            
        end       
    end
end

[x,y] = size(Ket60)
ii=0;
for CountC = 1:x
    for CountD = 1:length(Wavefunction)
        if Ket60(CountC,2:14) == Wavefunction(CountD,2:14)
            ii = ii +1;
            Ket60Val(ii,1) = Ket60(CountC,1) * Wavefunction(CountD,1)
        end
    end
end

[x,y] = size(Ket42)
ii=0;
for CountC = 1:x
    for CountD = 1:length(Wavefunction)
        if Ket42(CountC,2:14) == Wavefunction(CountD,2:14)
            ii = ii +1;
            Ket42Val(ii,1) = Ket42(CountC,1) * Wavefunction(CountD,1)
        end
    end
end

[x,y] = size(Ket24)
ii=0;
for CountC = 1:x
    for CountD = 1:length(Wavefunction)
        if Ket24(CountC,2:14) == Wavefunction(CountD,2:14)
            ii = ii +1;
            Ket24Val(ii,1) = Ket24(CountC,1) * Wavefunction(CountD,1)
        end
    end
end

[x,y] = size(Ket06)
ii=0;
for CountC = 1:x
    for CountD = 1:length(Wavefunction)
        if Ket06(CountC,2:14) == Wavefunction(CountD,2:14)
            ii = ii +1;
            Ket06Val(ii,1) = Ket06(CountC,1) * Wavefunction(CountD,1)
        end
    end
end
A = (sum(abs(Ket60Val.^2)))
B = (sum(abs(Ket42Val.^2)))
C = (sum(abs(Ket24Val.^2)))
D = (sum(abs(Ket06Val.^2)))
FidelityOverlap = A+B+C+D
% figure
% bar([A B C D])
% Vects
% toc
%% Creates Bat as in old method (LLL) plotting only most populated for each position

x = C0Valmax;
x1 = C2Valmax;
Eig = OrigCoeff;
%Eig = csvread('V_077721759.csv');
%%%%%%%%%%%%%%%%%%%%%%
Nt = N;
k=0;
LKetN0=zeros(1,10);
for N0 = 0
    for p = 0:N0
        k=k+1;
        
        Nk = nchoosek(N0,p);
        Norm = 1/((sqrt(factorial(N0))).*(sqrt(factorial(N-N0))));
        
        
        LKetN0(k,1) = Nk.*Norm.*(x^p).*(x1^(N0-p)).*(sqrt(factorial(p)))...
            .*sqrt(factorial(N0-p)).*sqrt(factorial(Nt-N0));
        LKetN0(k,2) = p;
        LKetN0(k,3) = 6-N0;
        LKetN0(k,4) = N0-p;
    end
end

k=0;
LKetN1 = zeros(2,10);
for N1 = 1
    for p = 0:N1
        k=k+1;

         Nk = nchoosek(N1,p)
         SumTerm = Nk.*((x)^(N1-p)).*((x1).^p)
         Norm = 1/((sqrt(factorial(N1))).*(sqrt(factorial(6-N1))))
         TotConst = SumTerm.*Norm.*Nk


        LKetN1(k,1) = Nk.*SumTerm.*Norm*(sqrt(factorial(p))*sqrt(factorial(6-N1))*sqrt(factorial(N1-p)));
        LKetN1(k,2) = p;
        LKetN1(k,3) = 6-N1;
        LKetN1(k,4) = N1-p;
    end
end

k=0;
LKetN2 = zeros(3,10);
for N2 = 2
    for p = 0:N2
        k=k+1;
        
        Nk = nchoosek(N2,p);
        Norm = 1/((sqrt(factorial(N2))).*(sqrt(factorial(Nt-N2))));
        
        
        LKetN2(k,1) = Nk.*Norm.*(x^p).*(x1^(N2-p)).*(sqrt(factorial(p)))...
            .*sqrt(factorial(N2-p)).*sqrt(factorial(Nt-N2));
        LKetN2(k,2) = p;
        LKetN2(k,3) = 6-N2;
        LKetN2(k,4) = N2-p;
    end
end


k=0;
LKetN3 = zeros(4,10);
for N3 = 3
    for p = 0:N3
        k=k+1;

         Nk = nchoosek(N3,p)
         SumTerm = Nk.*((x)^(N3-p)).*((x1).^p)
         Norm = 1/((sqrt(factorial(N3))).*(sqrt(factorial(6-N3))))
         TotConst = SumTerm.*Norm.*Nk


        LKetN3(k,1) = Nk.*SumTerm.*Norm*(sqrt(factorial(p))*sqrt(factorial(6-N3))*sqrt(factorial(N3-p)));
        LKetN3(k,2) = p;
        LKetN3(k,3) = 6-N3;
        LKetN3(k,4) = N3-p;
    end
end

k=0;
LKetN4 = zeros(5,10);
for N4 = 4
    for p = 0:N4
        k=k+1;
        
        Nk = nchoosek(N4,p);
        Norm = 1/((sqrt(factorial(N4))).*(sqrt(factorial(Nt-N4))));
        
        
        LKetN4(k,1) = Nk.*Norm.*(x^p).*(x1^(N4-p)).*(sqrt(factorial(p)))...
            .*sqrt(factorial(N4-p)).*sqrt(factorial(Nt-N4));
        LKetN4(k,2) = p;
        LKetN4(k,3) = 6-N4;
        LKetN4(k,4) = N4-p;
    end
end

k=0;
LKetN5 = zeros(6,10);
for N5 = 5
    for p = 0:N5
        k=k+1;

         Nk = nchoosek(N5,p)
         SumTerm = Nk.*((x)^(N5-p)).*((x1).^p)
         Norm = 1/((sqrt(factorial(N5))).*(sqrt(factorial(6-N5))))
         TotConst = SumTerm.*Norm.*Nk


        LKetN5(k,1) = Nk.*SumTerm.*Norm*(sqrt(factorial(p))*sqrt(factorial(6-N5))*sqrt(factorial(N5-p)));
        LKetN5(k,2) = p;
        LKetN5(k,3) = 6-N5;
        LKetN5(k,4) = N5-p;
    end
end


k=0;
LKetN6 = zeros(7,10);
for N6 = 6
    for p = 0:N6
        k=k+1;
        
        Nk = nchoosek(N6,p);
        Norm = 1/((sqrt(factorial(N6))).*(sqrt(factorial(Nt-N6))));
        
        
        LKetN6(k,1) = Nk.*Norm.*(x^p).*(x1^(N6-p)).*(sqrt(factorial(p)))...
            .*sqrt(factorial(N6-p)).*sqrt(factorial(Nt-N6));
        LKetN6(k,2) = p;
        LKetN6(k,3) = 6-N6;
        LKetN6(k,4) = N6-p;
    end
end

BasisChangeKet = [LKetN0;LKetN2;LKetN4;LKetN6];

%%%%%Above generates the basis transformation from |N_1,N_2> to the
%%%%%momentum basis used initially, next part finds fidelity with original
%%%%%wavefunction.


OrigStates = csvread('EvenBasisLL2.csv');
OrigWavefunction = [Eig(:,1) OrigStates];

%%%%%% Ket0 Part
l=0;m=0;
LKetN0SaveMat = zeros(10,1);
for CountA = 1:length(OrigWavefunction)
    l=l+1;
    if LKetN0(1,2:10) == OrigWavefunction(l,2:10)
        m=m+1;
        LKetN0SaveMat(m,1) = LKetN0(1,1)*OrigWavefunction(l,1) ;
    end
end

%%%%%%% Ket2 Part
l=0;m=0;
LKetN2SaveMat = zeros(10,1);
for CountA = 1:length(OrigWavefunction)
    l=l+1;
    if LKetN2(1,2:10) == OrigWavefunction(l,2:10)
        m=m+1;
        LKetN2SaveMat(m,1) = LKetN2(1,1)*OrigWavefunction(l,1);
    elseif LKetN2(2,2:10) == OrigWavefunction(l,2:10)
        m=m+1;
        LKetN2SaveMat(m,1) = LKetN2(2,1)*OrigWavefunction(l,1);
    elseif LKetN2(3,2:10) == OrigWavefunction(l,2:10)
        m=m+1;
        LKetN2SaveMat(m,1) = LKetN2(3,1)*OrigWavefunction(l,1);
    end
end

%%%%%%% Ket4 Part
l=0;m=0;
LKetN4SaveMat = zeros(10,1);
for CountA = 1:length(OrigWavefunction)
    l=l+1;
    if LKetN4(1,2:10) == OrigWavefunction(l,2:10)
        m=m+1;
        LKetN4SaveMat(m,1) = LKetN4(1,1)*OrigWavefunction(l,1);
    elseif LKetN4(2,2:10) == OrigWavefunction(l,2:10)
        m=m+1;
        LKetN4SaveMat(m,1) = LKetN4(2,1)*OrigWavefunction(l,1);
    elseif LKetN4(3,2:10) == OrigWavefunction(l,2:10)
        m=m+1;
        LKetN4SaveMat(m,1) = LKetN4(3,1)*OrigWavefunction(l,1);
    elseif LKetN4(4,2:10) == OrigWavefunction(l,2:10)
        m=m+1;
        LKetN4SaveMat(m,1) = LKetN4(4,1)*OrigWavefunction(l,1);
    elseif LKetN4(5,2:10) == OrigWavefunction(l,2:10)
        m=m+1;
        LKetN4SaveMat(m,1) = LKetN4(5,1)*OrigWavefunction(l,1);
    end
end

%%%%%Ket6 Part
l=0;m=0;
LKetN6SaveMat = zeros(10,1);
for CountA = 1:length(OrigWavefunction)
    l=l+1;
    if LKetN6(1,2:10) == OrigWavefunction(l,2:10)
        m=m+1;
        %KetN6(1,1);
        OrigWavefunction(l,1);
        LKetN6SaveMat(m,1) = LKetN6(1,1)*OrigWavefunction(l,1);
    elseif LKetN6(2,2:10) == OrigWavefunction(l,2:10)
        m=m+1;
        %KetN6(2,1);
        OrigWavefunction(l,1);
        LKetN6SaveMat(m,1) = LKetN6(2,1)*OrigWavefunction(l,1);
    elseif LKetN6(3,2:10) == OrigWavefunction(l,2:10)
        m=m+1;
        %KetN6(3,1);
        OrigWavefunction(l,1);
        LKetN6SaveMat(m,1) = LKetN6(3,1)*OrigWavefunction(l,1);
    elseif LKetN6(4,2:10) == OrigWavefunction(l,2:10)
        m=m+1;
        %KetN6(4,1);
        OrigWavefunction(l,1);
        LKetN6SaveMat(m,1) = LKetN6(4,1)*OrigWavefunction(l,1);
    elseif LKetN6(5,2:10) == OrigWavefunction(l,2:10)
        m=m+1;
        %KetN6(5,1);
        OrigWavefunction(l,1);
        LKetN6SaveMat(m,1) = LKetN6(5,1)*OrigWavefunction(l,1);
    elseif LKetN6(6,2:10) == OrigWavefunction(l,2:10)
        m=m+1;
        %KetN6(6,1);
        OrigWavefunction(l,1);
        LKetN6SaveMat(m,1) = LKetN6(6,1)*OrigWavefunction(l,1);
    elseif LKetN6(7,2:10) == OrigWavefunction(l,2:10)
        m=m+1;
        %KetN6(7,1);
        OrigWavefunction(l,1);
        LKetN6SaveMat(m,1) = LKetN6(7,1)*OrigWavefunction(l,1);
    end
end

LKetN0Overlap = (abs(sum(LKetN0SaveMat.^2)));
LKetN2Overlap = (abs(sum(LKetN2SaveMat.^2)));
LKetN4Overlap = (abs(sum(LKetN4SaveMat.^2)));
LKetN6Overlap = (abs(sum(LKetN6SaveMat.^2)));
OverlapsP = [LKetN0Overlap LKetN2Overlap LKetN4Overlap LKetN6Overlap];
figure
subplot(1,2,1)
bar([A B C D])
subplot(1,2,2)
bar([LKetN0Overlap LKetN2Overlap LKetN4Overlap LKetN6Overlap])




toc


















