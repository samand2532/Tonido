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

for Omega = 0.85;%0.7:0.001:0.9
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

[C0Valmax,C0Pos] = max(LargestVectTemp(:,2));
C0Valmax = LargestVectTemp(C0Pos,1) ;
LargestVectTemp(C0Pos,1) = 0;
LargestVectTemp(C0Pos,2) = 0;

[C2Valmax,C2Pos] = max(LargestVectTemp(:,2));
C2Valmax = LargestVectTemp(C2Pos,1)  ;
LargestVectTemp(C2Pos,1) = 0;
LargestVectTemp(C2Pos,2) = 0;

[C4Valmax,C4Pos] = max(LargestVectTemp(:,2));
C4Valmax = LargestVectTemp(C4Pos,1) ;
LargestVectTemp(C4Pos,1) = 0;
LargestVectTemp(C4Pos,2) = 0;

SecLargeVectTemp(:,1) = SecLargeVect(:);
SecLargeVectTemp(:,2) = abs(SecLargeVect(:));

[C1Valmax,C1Pos] = max(SecLargeVectTemp(:,2));
C1Valmax = SecLargeVectTemp(C1Pos,1)  ;
SecLargeVectTemp(C1Pos,1) = 0;
SecLargeVectTemp(C1Pos,2) = 0;

[C3Valmax,C3Pos] = max(SecLargeVectTemp(:,2));
C3Valmax = SecLargeVectTemp(C3Pos,1);
SecLargeVectTemp(C3Pos,1) = 0;
SecLargeVectTemp(C3Pos,2) = 0;

C2 = C0Valmax;
C4 = C2Valmax;
C6 = C4Valmax;
C1 = C1Valmax;
C3 = C3Valmax;

Psi1 = (C2^2)+(C4^2)+(C6^2);
Psi2 = (C1^2)+(C3^2);
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Arrangements of N1
for N1 = 0
    a = 0;
    for i = 0:6
        for j = 0:6
            for k = 0:6
                
                if i + j + k == N1
                    a = a + 1;
                    ijkMatN10(a,1) = N1;
                    ijkMatN10(a,2) = i;
                    ijkMatN10(a,3) = j;
                    ijkMatN10(a,4) = k;
                end
            end
        end
    end
end
for N1 = 1
    a = 0;
    for i = 0:6
        for j = 0:6
            for k = 0:6
                
                if (i + j + k) == N1
                    a = a + 1;
                    ijkMatN11(a,1) = N1;
                    ijkMatN11(a,2) = i;
                    ijkMatN11(a,3) = j;
                    ijkMatN11(a,4) = k;
                end
            end
        end
    end
end
for N1 = 2
    a = 0;
    for i = 0:6
        for j = 0:6
            for k = 0:6
                
                if i + j + k == N1
                    a = a + 1;
                    ijkMatN12(a,1) = N1;
                    ijkMatN12(a,2) = i;
                    ijkMatN12(a,3) = j;
                    ijkMatN12(a,4) = k;
                end
            end
        end
    end
end
for N1 = 3
    a = 0;
    for i = 0:6
        for j = 0:6
            for k = 0:6
                
                if i + j + k == N1
                    a = a + 1;
                    ijkMatN13(a,1) = N1;
                    ijkMatN13(a,2) = i;
                    ijkMatN13(a,3) = j;
                    ijkMatN13(a,4) = k;
                end
            end
        end
    end
end
for N1 = 4
    a = 0;
    for i = 0:6
        for j = 0:6
            for k = 0:6
                
                if i + j + k == N1
                    a = a + 1;
                    ijkMatN14(a,1) = N1;
                    ijkMatN14(a,2) = i;
                    ijkMatN14(a,3) = j;
                    ijkMatN14(a,4) = k;
                end
            end
        end
    end
end
for N1 = 5
    a = 0;
    for i = 0:6
        for j = 0:6
            for k = 0:6
                
                if i + j + k == N1
                    a = a + 1;
                    ijkMatN15(a,1) = N1;
                    ijkMatN15(a,2) = i;
                    ijkMatN15(a,3) = j;
                    ijkMatN15(a,4) = k;
                end
            end
        end
    end
end
for N1 = 6
    a = 0;
    for i = 0:6
        for j = 0:6
            for k = 0:6
                
                if i + j + k == N1
                    a = a + 1;
                    ijkMatN16(a,1) = N1;
                    ijkMatN16(a,2) = i;
                    ijkMatN16(a,3) = j;
                    ijkMatN16(a,4) = k;
                end
            end
        end
    end
end
ijkOut = [ijkMatN10; ijkMatN11; ijkMatN12; ijkMatN13; ijkMatN14; ijkMatN15; ijkMatN16];
bb = 0;
for CountA = 1:length(ijkOut)
    N1 = ijkOut(CountA,1);
    N2 = (N - N1);
    i = ijkOut(CountA,2);
    j = ijkOut(CountA,3);
    k = ijkOut(CountA,4);
    for CountB = N2
        for p = 0:N2
            q = N2 - p;
            bb = bb + 1;
            BiTriOut(bb,1)=N1;
            BiTriOut(bb,2)=N2;
            BiTriOut(bb,C0Pos+3)=i;
            BiTriOut(bb,C2Pos+3)=j;
            BiTriOut(bb,C4Pos+3)=k;
            BiTriOut(bb,C1Pos+3)=p;
            BiTriOut(bb,C3Pos+3)=q;
            
        end
    end
end
BiTriVect = zeros(1,16);
cc=0;
%%%% Section stops Mt = -1 and mt of n=1 occuring more than once.
for BiTriVectCount = 1:length(BiTriOut)
    
    N1 = BiTriOut(BiTriVectCount,1);
    N2 = BiTriOut(BiTriVectCount,2);
    i = BiTriOut(BiTriVectCount,C0Pos+3);
    j = BiTriOut(BiTriVectCount,C2Pos+3);
    k = BiTriOut(BiTriVectCount,C4Pos+3);
    p = BiTriOut(BiTriVectCount,C1Pos+3);
    q = BiTriOut(BiTriVectCount,C3Pos+3);
    if (p < 2) && (k < 2)
        cc = cc + 1;
        BiTriVect(cc,1) = N1;
        BiTriVect(cc,2) = N2;
        BiTriVect(cc,C1Pos+3) = p;
        BiTriVect(cc,C0Pos+3) = i;
        BiTriVect(cc,C3Pos+3) = q;
        BiTriVect(cc,C2Pos+3) = j;
        BiTriVect(cc,C4Pos+3) = k;
    end
end
for ConstbiTri = 1:length(BiTriVect)
    N1 = BiTriVect(ConstbiTri,1);
    N2 = BiTriVect(ConstbiTri,2);
    i = BiTriVect(ConstbiTri,C0Pos+3);
    j = BiTriVect(ConstbiTri,C2Pos+3);
    k = BiTriVect(ConstbiTri,C4Pos+3);
    p = BiTriVect(ConstbiTri,C1Pos+3);
    q = BiTriVect(ConstbiTri,C3Pos+3);
    
    BiTriVect(ConstbiTri,3) = (1/((sqrt(factorial(N1)))*(sqrt(factorial(N2)))))...
        *(factorial(N1)/(factorial(i)*factorial(j)*factorial(k)))...
        *(C2^i)*(C4^j)*(C6^k)*sqrt(factorial(i))*sqrt(factorial(j))*sqrt(factorial(k))...
        *nchoosek(N2,p)*(C1^p)*(C3^q)*sqrt(factorial(p))*sqrt(factorial(q));
end

Wavefunction;
cc = 0;

KetN60 = 0;
KetN51 = 0;
KetN42 = 0;
KetN33 = 0;
KetN24 = 0;
KetN15 = 0;
KetN06 = 0;

for BiTriVectCount = 1:length(BiTriVect)
    
    if BiTriVect(BiTriVectCount,1) == 6
        
        for CountWave = 1:length(Wavefunction)
            if BiTriVect(BiTriVectCount,4:16) == Wavefunction(CountWave,2:14)
                cc = cc + 1;
                KetN60(cc,1) = BiTriVect(BiTriVectCount,3)*Wavefunction(CountWave,1);
           
            end
        end
        
        
    elseif BiTriVect(BiTriVectCount,1) == 5
        for CountWave = 1:length(Wavefunction)
            if BiTriVect(BiTriVectCount,4:16) == Wavefunction(CountWave,2:14)
                cc = cc + 1;
                KetN51(cc,1) = BiTriVect(BiTriVectCount,3)*Wavefunction(CountWave,1);
            
            end
        end
        
    elseif BiTriVect(BiTriVectCount,1) == 4
        for CountWave = 1:length(Wavefunction)
            if BiTriVect(BiTriVectCount,4:16) == Wavefunction(CountWave,2:14)
                cc = cc + 1;
                KetN42(cc,1) = BiTriVect(BiTriVectCount,3)*Wavefunction(CountWave,1);
            
            end
        end
        
        
    elseif BiTriVect(BiTriVectCount,1) == 3
        for CountWave = 1:length(Wavefunction)
            if BiTriVect(BiTriVectCount,4:16) == Wavefunction(CountWave,2:14)
                cc = cc + 1;
                KetN33(cc,1) = BiTriVect(BiTriVectCount,3)*Wavefunction(CountWave,1);
             
            end
        end
        
        
    elseif BiTriVect(BiTriVectCount,1) == 2
        for CountWave = 1:length(Wavefunction)
            if BiTriVect(BiTriVectCount,4:16) == Wavefunction(CountWave,2:14)
                cc = cc + 1;
                KetN24(cc,1) = BiTriVect(BiTriVectCount,3)*Wavefunction(CountWave,1);
            
            end
        end
        
        
    elseif BiTriVect(BiTriVectCount,1) == 1
        for CountWave = 1:length(Wavefunction)
            if BiTriVect(BiTriVectCount,4:16) == Wavefunction(CountWave,2:14)
                cc = cc + 1;
                KetN15(cc,1) = BiTriVect(BiTriVectCount,3)*Wavefunction(CountWave,1);
            
            end
        end
        
        
        
    elseif BiTriVect(BiTriVectCount,1) == 0
        for CountWave = 1:length(Wavefunction)
            if BiTriVect(BiTriVectCount,4:16) == Wavefunction(CountWave,2:14)
                cc = cc + 1;
                KetN06(cc,1) = BiTriVect(BiTriVectCount,3)*Wavefunction(CountWave,1);
            
            end
        end
        
        
    end
end
Overlap60 = sum(KetN60).^2;
Overlap51 = sum(KetN51).^2;
Overlap42 = sum(KetN42).^2;
Overlap33 = sum(KetN33).^2;
Overlap24 = sum(KetN24).^2;
Overlap15 = sum(KetN15).^2;
Overlap06 = sum(KetN06).^2;




Overlaps = [ Overlap60 Overlap51 Overlap42 Overlap33 Overlap24 Overlap15 Overlap06 ]
bar(Overlaps)
toc