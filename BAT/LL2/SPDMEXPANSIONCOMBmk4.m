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
format Short
%%%%%%%%%% What fidelity we want to achieve
Tolerance = 0.99
%%% T = 0.75 gives C0 and C1
%%% T = 0.85 gives C0 C1 and C3 
%%% T = 0.97 gives C0 C2 C1 and C3
%%%%%%%%%%
for Omega = 0.82856708053820;%0.001:0.9
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
%%%%% if largest coeff doesnt exceed tolerence find the next coefficiant

[C2Valmax,C2Pos] = max(LargestVectTemp(:,2));
C2Valmax = LargestVectTemp(C2Pos,1);
LargestVectTemp(C2Pos,1) = 0;
LargestVectTemp(C2Pos,2) = 0;

[C4Valmax,C4Pos] = max(LargestVectTemp(:,2));
C4Valmax = LargestVectTemp(C4Pos,1);
LargestVectTemp(C4Pos,1) = 0;
LargestVectTemp(C4Pos,2) = 0;


TempSecLargeVect(:,1) = SecLargeVect(:);
TempSecLargeVect(:,2) = abs(SecLargeVect(:));

[C1Valmax,C1Pos] = max(TempSecLargeVect(:,2));
C1Valmax = TempSecLargeVect(C1Pos,1);
TempSecLargeVect(C1Pos,1) = 0;
TempSecLargeVect(C1Pos,2) = 0;

[C3Valmax,C3Pos] = max(TempSecLargeVect(:,2));
C3Valmax = TempSecLargeVect(C3Pos,1);
TempSecLargeVect(C3Pos,1) = 0;
TempSecLargeVect(C3Pos,2) = 0;
% disp('C0')
% [C0Valmax C0Pos]
% disp('C2')
% [C2Valmax C2Pos]
% disp('C4')
% [C4Valmax C4Pos]
% disp('C1')
% [C1Valmax C1Pos]
% disp('C3')
% [C3Valmax C3Pos]

if C0Valmax.^2 > Tolerance
    C0Coeff = C0Valmax;
    C0Pos = C0Pos;
elseif (C0Valmax.^2)+(C2Valmax.^2)
    C0Coeff = C0Valmax
    C0Pos = C0Pos;
    C2Coeff = C2Valmax;
    C2Pos = C2Pos;
else C0Coeff = C0Valmax
    C0Pos = C0Pos;
    C2Coeff = C2Valmax;
    C2Pos = C2Pos;
    C4Coeff = C4Valmax;
    C4Pos = C4Pos;
end
if C1Valmax.^2 > Tolerance
    C1Coeff = C1Valmax;
    C1Pos = C1Pos;
else 
    C1Coeff = C1Valmax;
    C1Pos = C1Pos;
    C3Coeff = C3Valmax;
    C3Pos = C3Pos;
end

%%% Exist command return 1 if the variable tested exists
C0Exist = exist('C0Coeff');
C2Exist = exist('C2Coeff');
C4Exist = exist('C4Coeff');
C1Exist = exist('C1Coeff');
C3Exist = exist('C3Coeff');

[LargestVect SecLargeVect];
%Finds Psi, but only if the other Coefficiants do not exist
if C0Exist == 1 && C2Exist == 0 && C4Exist == 0
    Psi1 = (C0Coeff^2)
    disp('C0')
    C0Pos
elseif C0Exist == 1 && C2Exist == 1 && C4Exist == 0
    Psi1 = (C0Coeff^2)+(C2Coeff^2)
    disp('C0')
    disp('C2')
    C0Pos
    C2Pos
elseif C0Exist == 1 && C2Exist == 1 && C4Exist == 1
    disp('C0')
    disp('C2')
    disp('C4')
    C0Pos
    C2Pos
    C4Pos
    Psi1 = (C0Coeff^2)+(C2Coeff^2)+(C4Coeff^2)
end

if C1Exist == 1 && C3Exist == 0
    Psi2 = (C1Coeff^2)
    disp('C1')
    C1Pos
else Psi2 = (C1Coeff^2) + (C3Coeff^2)
    disp('C1')
    disp('C3')
    C1Pos
    C3Pos
end

%%%%%%%%%
Ket06 = zeros(1,23);
Ket15 = zeros(1,23);
Ket24 = zeros(1,23);
Ket33 = zeros(1,23);
Ket42 = zeros(1,23);
Ket51 = zeros(1,23);
Ket60 = zeros(1,23);
%% C0 and C1  exists 
if C0Exist==1 && C2Exist==0 && C4Exist ==0 && C1Exist==1 && C3Exist==0
    %%%singlexsingle
    
    N1 = 0;
    N2 = 6;
    
    Ket06(1,1) = 1/((sqrt(factorial(N1)))*(sqrt(factorial(N2))))...
        *(C0Coeff^N1)*(C1Coeff^N2)*sqrt(factorial(N1))*sqrt(factorial(N2));
    Ket06(1,C0Pos+1) = N1;
    Ket06(1,C1Pos+1) = N2;
    
    N1 = 2;
    N2 = 4;
    
    Ket24(1,1) = 1/((sqrt(factorial(N1)))*(sqrt(factorial(N2))))...
        *(C0Coeff^N1)*(C1Coeff^N2)*sqrt(factorial(N1))*sqrt(factorial(N2));
    Ket24(1,C0Pos+1) = N1;
    Ket24(1,C1Pos+1) = N2;
    
    
    N1 = 4;
    N2 = 2;
    
    Ket42(1,1) = 1/((sqrt(factorial(N1)))*(sqrt(factorial(N2))))...
        *(C0Coeff^N1)*(C1Coeff^N2)*sqrt(factorial(N1))*sqrt(factorial(N2));
    Ket42(1,C0Pos+1) = N1;
    Ket42(1,C1Pos+1) = N2;
    
    
    N1 = 6;
    N2 = 0;
    
    Ket60(1,1) = 1/((sqrt(factorial(N1)))*(sqrt(factorial(N2))))...
        *(C0Coeff^N1)*(C1Coeff^N2)*sqrt(factorial(N1))*sqrt(factorial(N2));
    
    Ket60(1,C0Pos+1) = N1;
    Ket60(1,C1Pos+1) = N2;
    

%% C0,C1 and C3 exists    
elseif C0Exist==1 && C2Exist==0 && C4Exist ==0 && C1Exist==1 && C3Exist==1
   
    for N1 = 0;
        N2 = N - N1;
        k = 0;
        for p = 0:N2
            k = k + 1;
            Ket06(k,1) = 1/((sqrt(factorial(N1)))*(sqrt(factorial(N2))))...
                *nchoosek(N2,p)...
                *(C0Coeff^N1)*(C1Coeff^p)*(C0Coeff^(N2-p))...
                *sqrt(factorial(N1))*sqrt(factorial(p))*sqrt(factorial(N2-p));
            
            Ket06(k,C0Pos+1) = N1;
            Ket06(k,C1Pos+1) = p;
            Ket06(k,C3Pos+1) = (N2-p);
        end
        
        
    end
    
    for N1 = 1;
        N2 = N - N1;
        k = 0;
        for p = 0:N2
            k = k + 1;
            Ket15(k,1) = 1/(sqrt(factorial(N1)))*(sqrt(factorial(N2)))...
                *nchoosek(N2,p)...
                *(C0Coeff^N1)*(C1Coeff^p)*(C0Coeff^(N2-p))...
                *sqrt(factorial(N1))*sqrt(factorial(p))*sqrt(factorial(N2-p));
            
            Ket15(k,C0Pos+1) = N1;
            Ket15(k,C1Pos+1) = p;
            Ket15(k,C3Pos+1) = (N2-p);
        end
    end
    
    for N1 = 2;
        N2 = N - N1;
        k = 0;
        for p = 0:N2
            k = k + 1;
            Ket24(k,1) = 1/(sqrt(factorial(N1)))*(sqrt(factorial(N2)))...
                *nchoosek(N2,p)...
                *(C0Coeff^N1)*(C1Coeff^p)*(C0Coeff^(N2-p))...
                *sqrt(factorial(N1))*sqrt(factorial(p))*sqrt(factorial(N2-p));
            
            Ket24(k,C0Pos+1) = N1;
            Ket24(k,C1Pos+1) = p;
            Ket24(k,C3Pos+1) = (N2-p);
        end
    end
    
    for N1 = 3;
        N2 = N - N1;
        k = 0;
        for p = 0:N2
            k = k + 1;
            Ket33(k,1) = 1/(sqrt(factorial(N1)))*(sqrt(factorial(N2)))...
                *nchoosek(N2,p)...
                *(C0Coeff^N1)*(C1Coeff^p)*(C0Coeff^(N2-p))...
                *sqrt(factorial(N1))*sqrt(factorial(p))*sqrt(factorial(N2-p));
            
            Ket33(k,C0Pos+1) = N1;
            Ket33(k,C1Pos+1) = p;
            Ket33(k,C3Pos+1) = (N2-p);
        end
    end
    
    for N1 = 4;
        N2 = N - N1;
        k = 0;
        for p = 0:N2
            k = k + 1;
            Ket42(k,1) = 1/(sqrt(factorial(N1)))*(sqrt(factorial(N2)))...
                *nchoosek(N2,p)...
                *(C0Coeff^N1)*(C1Coeff^p)*(C0Coeff^(N2-p))...
                *sqrt(factorial(N1))*sqrt(factorial(p))*sqrt(factorial(N2-p));
            
            Ket42(k,C0Pos+1) = N1;
            Ket42(k,C1Pos+1) = p;
            Ket42(k,C3Pos+1) = (N2-p);
        end
    end
    
    for N1 = 5;
        N2 = N - N1;
        k = 0;
        for p = 0:N2
            k = k + 1;
            Ket51(k,1) = 1/(sqrt(factorial(N1)))*(sqrt(factorial(N2)))...
                *nchoosek(N2,p)...
                *(C0Coeff^N1)*(C1Coeff^p)*(C0Coeff^(N2-p))...
                *sqrt(factorial(N1))*sqrt(factorial(p))*sqrt(factorial(N2-p));
            
            Ket51(k,C0Pos+1) = N1;
            Ket51(k,C1Pos+1) = p;
            Ket51(k,C3Pos+1) = (N2-p);
        end
    end
    
    for N1 = 6;
        N2 = N - N1;
        k = 0;
        for p = 0:N2
            k = k + 1;
            Ket60(k,1) = 1/(sqrt(factorial(N1)))*(sqrt(factorial(N2)))...
                *nchoosek(N2,p)...
                *(C0Coeff^N1)*(C1Coeff^p)*(C0Coeff^(N2-p))...
                *sqrt(factorial(N1))*sqrt(factorial(p))*sqrt(factorial(N2-p));
            
            Ket60(k,C0Pos+1) = N1;
            Ket60(k,C1Pos+1) = p;
            Ket60(k,C3Pos+1) = (N2-p);
        end
    end
    
    
    
  %% End Single Bi
  %% Start of BixSingle
elseif C0Exist==1 && C2Exist==1 && C4Exist ==0 && C1Exist==1 && C3Exist==0
     %%%bixsingle
     
     
    
    %% end of BixSingle  
    %% BixBi
elseif C0Exist==1 && C2Exist==1 && C4Exist ==0 && C1Exist==1 && C3Exist==1
    %%%bixbi
    for N1 = 0
        N2 = N - N1;
        %%%% Starting with coeff from largest first.
        k = 0;
        for p = 0:N1
            for q = 0:N2
                k = k +1;
                Ket06(k,1) = 1/(sqrt(factorial(N1)))*(sqrt(factorial(N2)))...
                    *nchoosek(N1,p)*nchoosek(N2,q)...
                    *(C0Coeff^p)*(C2Coeff^(N1-p))*(C1Coeff^q)*(C3Coeff^(N2-q))...
                    *sqrt(factorial(p))*sqrt(factorial(N1-p))...
                    *sqrt(factorial(q))*sqrt(factorial(N2-q));
                nchoosek(N1,p)
                sqrt(factorial(p))
                Ket06(k,C0Pos+1) = p;
                Ket06(k,C2Pos+1) = (N1-p);
                Ket06(k,C1Pos+1) = q;
                Ket06(k,C3Pos+1) = (N2-q);
            end
        end
        
        
        
    end
    
    for N1 = 1
        N2 = N - N1;
        k = 0;
        for p = 0:N1
            for q = 0:N2
                k = k +1;
                Ket15(k,1) = 1/(sqrt(factorial(N1)))*(sqrt(factorial(N2)))...
                    *nchoosek(N1,p)*nchoosek(N2,q)...
                    *(C0Coeff^p)*(C2Coeff^(N1-p))*(C1Coeff^q)*(C3Coeff^(N2-q))...
                    *sqrt(factorial(p))*sqrt(factorial(N1-p))...
                    *sqrt(factorial(q))*sqrt(factorial(N2-q));
                
                Ket15(k,C0Pos+1) = p;
                Ket15(k,C2Pos+1) = (N1-p);
                Ket15(k,C1Pos+1) = q;
                Ket15(k,C3Pos+1) = (N2-q);
            end
        end
    end
    
    for N1 = 2
        N2 = N - N1;
        k = 0;
        for p = 0:N1
            for q = 0:N2
                k = k +1;
                Ket24(k,1) = 1/(sqrt(factorial(N1)))*(sqrt(factorial(N2)))...
                    *nchoosek(N1,p)*nchoosek(N2,q)...
                    *(C0Coeff^p)*(C2Coeff^(N1-p))*(C1Coeff^q)*(C3Coeff^(N2-q))...
                    *sqrt(factorial(p))*sqrt(factorial(N1-p))...
                    *sqrt(factorial(q))*sqrt(factorial(N2-q));
                
                Ket24(k,C0Pos+1) = p;
                Ket24(k,C2Pos+1) = (N1-p);
                Ket24(k,C1Pos+1) = q;
                Ket24(k,C3Pos+1) = (N2-q);
            end
        end
    end
    
    for N1 = 3
        N2 = N - N1;
        k = 0;
        for p = 0:N1
            for q = 0:N2
                k = k +1;
                Ket33(k,1) = 1/(sqrt(factorial(N1)))*(sqrt(factorial(N2)))...
                    *nchoosek(N1,p)*nchoosek(N2,q)...
                    *(C0Coeff^p)*(C2Coeff^(N1-p))*(C1Coeff^q)*(C3Coeff^(N2-q))...
                    *sqrt(factorial(p))*sqrt(factorial(N1-p))...
                    *sqrt(factorial(q))*sqrt(factorial(N2-q));
                
                Ket33(k,C0Pos+1) = p;
                Ket33(k,C2Pos+1) = (N1-p);
                Ket33(k,C1Pos+1) = q;
                Ket33(k,C3Pos+1) = (N2-q);
            end
        end
    end
    
    for N1 = 4
        N2 = N - N1;
        k = 0;
        for p = 0:N1
            for q = 0:N2
                k = k +1;
                Ket42(k,1) = 1/(sqrt(factorial(N1)))*(sqrt(factorial(N2)))...
                    *nchoosek(N1,p)*nchoosek(N2,q)...
                    *(C0Coeff^p)*(C2Coeff^(N1-p))*(C1Coeff^q)*(C3Coeff^(N2-q))...
                    *sqrt(factorial(p))*sqrt(factorial(N1-p))...
                    *sqrt(factorial(q))*sqrt(factorial(N2-q));
                
                Ket42(k,C0Pos+1) = p;
                Ket42(k,C2Pos+1) = (N1-p);
                Ket42(k,C1Pos+1) = q;
                Ket42(k,C3Pos+1) = (N2-q);
            end
        end
    end
    
    for N1 = 5
        N2 = N - N1;
        k = 0;
        for p = 0:N1
            for q = 0:N2
                k = k +1;
                Ket51(k,1) = 1/(sqrt(factorial(N1)))*(sqrt(factorial(N2)))...
                    *nchoosek(N1,p)*nchoosek(N2,q)...
                    *(C0Coeff^p)*(C2Coeff^(N1-p))*(C1Coeff^q)*(C3Coeff^(N2-q))...
                    *sqrt(factorial(p))*sqrt(factorial(N1-p))...
                    *sqrt(factorial(q))*sqrt(factorial(N2-q));
                
                Ket51(k,C0Pos+1) = p;
                Ket51(k,C2Pos+1) = (N1-p);
                Ket51(k,C1Pos+1) = q;
                Ket51(k,C3Pos+1) = (N2-q);
            end
        end
    end
    
    for N1 = 6
        N2 = N - N1;
        k = 0;
        for p = 0:N1
            for q = 0:N2
                k = k +1;
                Ket60(k,1) = 1/(sqrt(factorial(N1)))*(sqrt(factorial(N2)))...
                    *nchoosek(N1,p)*nchoosek(N2,q)...
                    *(C0Coeff^p)*(C2Coeff^(N1-p))*(C1Coeff^q)*(C3Coeff^(N2-q))...
                    *sqrt(factorial(p))*sqrt(factorial(N1-p))...
                    *sqrt(factorial(q))*sqrt(factorial(N2-q));
                
                Ket60(k,C0Pos+1) = p;
                Ket60(k,C2Pos+1) = (N1-p);
                Ket60(k,C1Pos+1) = q;
                Ket60(k,C3Pos+1) = (N2-q);
            end
        end
    end
    
    %% end BixBi
elseif C0Exist==1 && C2Exist==1 && C4Exist ==1 && C1Exist==1 && C3Exist==0
    %%%trixsingle
elseif C0Exist==1 && C2Exist==1 && C4Exist ==1 && C1Exist==1 && C3Exist==1
    %%%trixbi
end

%% finds cross of TM with Wavefunction
[x,y] = size(Ket60);
for CountA = 1:x
    for CountB = 1:length(Wavefunction)
        if Ket60(CountA,2:10) == Wavefunction(CountB,2:10)
            Ket60;
            Ket60Val(1,1) = Ket60(CountA,1)*Wavefunction(CountB,1);
        end
    end
end

[x,y] = size(Ket42);
for CountA = 1:x
    for CountB = 1:length(Wavefunction)
        if Ket42(CountA,2:10) == Wavefunction(CountB,2:10)
            Ket42;
            Ket42Val(1,1) = Ket42(CountA,1)*Wavefunction(CountB,1);
        end
    end
end

[x,y] = size(Ket24);
for CountA = 1:x
    for CountB = 1:length(Wavefunction)
        if Ket24(CountA,2:10) == Wavefunction(CountB,2:10)
            Ket24;
            Ket24Val(1,1) = Ket24(CountA,1)*Wavefunction(CountB,1);
        end
    end
end


[x,y] = size(Ket06);
for CountA = 1:x
    for CountB = 1:length(Wavefunction)
        if Ket06(CountA,2:10) == Wavefunction(CountB,2:10)
            Ket06;
            Ket06Val(1,1) = Ket06(CountA,1)*Wavefunction(CountB,1);
        end
    end
end
[LargestVect SecLargeVect]
Overlaps = [Ket60Val.^2 Ket24Val.^2 Ket42Val.^2 Ket06Val.^2]
Fidelity = sum(Overlaps)
bar(Overlaps)
toc