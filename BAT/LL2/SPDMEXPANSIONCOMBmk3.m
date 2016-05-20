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
Tolerance = 0.6
%%%%%%%%%%
for Omega = 0.85;%0.7:0.001:0.9
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
%[ LargestVect  SecLargeVect]
%%%% finds the 3 most populated states of the largestVect

LargestVectTemp(:,1) = LargestVect(:);
LargestVectTemp(:,2) = abs(LargestVect(:));

[C0Valmax,C0Pos] = max(LargestVectTemp(:,2));
C0Valmax = LargestVectTemp(C0Pos,1) ;
LargestVectTemp(C0Pos,1) = 0;
LargestVectTemp(C0Pos,2) = 0;
%%%%% if largest coeff doesnt exceed tolerence find the next coefficiant
if C0Valmax^2 > Tolerance
    C0Coeff = C0Valmax;
    C0Pos = C0Pos;
else
    C0Coeff = C0Valmax;
    C0Pos = C0Pos;
    [C2Valmax,C2Pos] = max(LargestVectTemp(:,2));
    C2Valmax = LargestVectTemp(C2Pos,1)  ;
    LargestVectTemp(C2Pos,1) = 0;
    LargestVectTemp(C2Pos,2) = 0;
    
    if (C0Valmax^2 + C2Valmax^2) > Tolerance
        C2Coeff = C2Valmax;
        C2Pos = C2Pos;
    else
        [C4Valmax,C4Pos] = max(LargestVectTemp(:,2));
        C4Valmax = LargestVectTemp(C4Pos,1) ;
        LargestVectTemp(C4Pos,1) = 0;
        LargestVectTemp(C4Pos,2) = 0;
        C4Coeff = C4Valmax;
        C4Pos = C4Pos;
    end
end
SecLargeVectTemp(:,1) = SecLargeVect(:);
SecLargeVectTemp(:,2) = abs(SecLargeVect(:));

[C1Valmax,C1Pos] = max(SecLargeVectTemp(:,2));
C1Valmax = SecLargeVectTemp(C1Pos,1)  ;
SecLargeVectTemp(C1Pos,1) = 0;
SecLargeVectTemp(C1Pos,2) = 0;

if (C1Valmax^2) > Tolerance
    C1Coeff = C1Valmax;
    C1Pos = C1Pos;
else
    C1Coeff = C1Valmax;
    C1Pos = C1Pos;
    [C3Valmax,C3Pos] = max(SecLargeVectTemp(:,2));
    C3Valmax = SecLargeVectTemp(C3Pos,1);
    SecLargeVectTemp(C3Pos,1) = 0;
    SecLargeVectTemp(C3Pos,2) = 0;
    C3Coeff = C3Valmax;
    C3Pos = C3Pos;
end
%%% Exist command return 1 if the variable tested exists
C0Exist = exist('C0Coeff');
C2Exist = exist('C2Coeff');
C4Exist = exist('C4Coeff');
C1Exist = exist('C1Coeff');
C3Exist = exist('C3Coeff');

[LargestVect SecLargeVect]
%Finds Psi, but only if the other Coefficiants do not exist
if C0Exist == 1 && C2Exist == 0 && C4Exist == 0
    Psi1 = (C0Coeff^2)
    disp('C0')
elseif C0Exist == 1 && C2Exist == 1 && C4Exist == 0
    Psi1 = (C0Coeff^2)+(C2Coeff^2)
    disp('C0')
    disp('C2')
elseif C0Exist == 1 && C2Exist == 1 && C4Exist == 1
    disp('C0')
    disp('C2')
    disp('C4')
    Psi1 = (C0Coeff^2)+(C2Coeff^2)+(C4Coeff^2)
end

if C1Exist == 1 && C3Exist == 0
    Psi2 = (C1Coeff^2)
    disp('C1')
else Psi2 = (C1Coeff^2) + (C3Coeff^2)
    disp('C1')
    disp('C3')
end

%%%%%%%%%
Ket06 = zeros(1,25);
Ket15 = zeros(1,25);
Ket24 = zeros(1,25);
Ket33 = zeros(1,25);
Ket42 = zeros(1,25);
Ket51 = zeros(1,25);
Ket60 = zeros(1,25);
%% C0 and C1  exists 
if C0Exist==1 && C2Exist==0 && C4Exist ==0 && C1Exist==1 && C3Exist==0
    %%%singlexsingle
    k=1;
    N1 = 0;
    N2 = 6;
    
    Ket06(k,3) = 1/(sqrt(factorial(N1))*sqrt(factorial(N2)))...
    *(C0Coeff^N1)*(C1Coeff^N2)*sqrt(factorial(N1))*sqrt(factorial(N2));
    Ket06(k,C0Pos+3) = N1;
    Ket06(k,C1Pos+3) = N2;
    
    N1 = 2;
    N2 = 4;
    
    Ket24(k,3) = 1/(sqrt(factorial(N1))*sqrt(factorial(N2)))...
    *(C0Coeff^N1)*(C1Coeff^N2)*sqrt(factorial(N1))*sqrt(factorial(N2));
    Ket24(k,C0Pos+3) = N1;
    Ket24(k,C1Pos+3) = N2;
    
    N1 = 4;
    N2 = 2;
    
    Ket42(k,3) = 1/(sqrt(factorial(N1))*sqrt(factorial(N2)))...
    *(C0Coeff^N1)*(C1Coeff^N2)*sqrt(factorial(N1))*sqrt(factorial(N2));
    Ket42(k,C0Pos+3) = N1;
    Ket42(k,C1Pos+3) = N2;
    
    N1 = 6;
    N2 = 0;
    
    Ket60(k,3) = 1/(sqrt(factorial(N1))*sqrt(factorial(N2)))...
    *(C0Coeff^N1)*(C1Coeff^N2)*sqrt(factorial(N1))*sqrt(factorial(N2));
    Ket60(k,C0Pos+3) = N1;
    Ket60(k,C1Pos+3) = N2;
    %%%%% finds overlaps with wavefunction
    
    
    
    
%% C0,C1 and C3 exists    
elseif C0Exist==1 && C2Exist==0 && C4Exist ==0 && C1Exist==1 && C3Exist==1
    %% singlexbi
    for N1 = 0:6
        N1;
        N2 = 6 - N1;
        [N1 N2];
        
        if N1 == 0;
            k=0;
            for p = 0:N2
                N1;
                N2;
                k=k+1;
                Ket06(k,3) = 1/(sqrt(factorial(N1))*sqrt(factorial(N2)))...
                    *nchoosek(N2,p)*(C0Coeff^N1)*(C1Coeff^(N2-p))*(C3Coeff^p)...
                    *sqrt(factorial(N1))*sqrt(factorial(N2-p))*sqrt(factorial(p));
                Ket06(k,1) = N1;
                Ket06(k,2) = N2;
                %KetN1_0(k,C0Pos+3) = N1;
                Ket06(k,C1Pos+4) = N2-p;
                Ket06(k,C3Pos+4) = p;
                
            end
        elseif N1 == 1;
            k=0;
            for p = 0:N2
                N1;
                N2;
                k=k+1;
                Ket15(k,3) = 1/(sqrt(factorial(N1))*sqrt(factorial(N2)))...
                    *nchoosek(N2,p)*(C0Coeff^N1)*(C1Coeff^(N2-p))*(C3Coeff^p)...
                    *sqrt(factorial(N1))*sqrt(factorial(N2-p))*sqrt(factorial(p));
                Ket15(k,1) = N1;
                Ket15(k,2) = N2;
                %KetN1_1(k,C0Pos+3) = N1;
                Ket15(k,C1Pos+4) = N2-p;
                Ket15(k,C3Pos+4) = p;
            end
        elseif N1 == 2;
            k=0;
            for p = 0:N2
                N1;
                N2;
                k=k+1;
                Ket24(k,3) = 1/(sqrt(factorial(N1))*sqrt(factorial(N2)))...
                    *nchoosek(N2,p)*(C0Coeff^N1)*(C1Coeff^(N2-p))*(C3Coeff^p)...
                    *sqrt(factorial(N1))*sqrt(factorial(N2-p))*sqrt(factorial(p));
                Ket24(k,1) = N1;
                Ket24(k,2) = N2;
                %KetN1_2(k,C0Pos+3) = N1;
                Ket24(k,C1Pos+4) = N2-p;
                Ket24(k,C3Pos+4) = p;
            end
        elseif N1 == 3;
            k=0;
            for p = 0:N2
                N1;
                N2;
                k=k+1;
                Ket33(k,3) = 1/(sqrt(factorial(N1))*sqrt(factorial(N2)))...
                    *nchoosek(N2,p)*(C0Coeff^N1)*(C1Coeff^(N2-p))*(C3Coeff^p)...
                    *sqrt(factorial(N1))*sqrt(factorial(N2-p))*sqrt(factorial(p));
                Ket33(k,1) = N1;
                Ket33(k,2) = N2;
                %KetN1_3(k,C0Pos+3) = N1;
                Ket33(k,C1Pos+4) = N2-p;
                Ket33(k,C3Pos+4) = p;
            end
            
        elseif N1 == 4;
            k=0;
            for p = 0:N2
                N1;
                N2;
                k=k+1;
                Ket42(k,3) = 1/(sqrt(factorial(N1))*sqrt(factorial(N2)))...
                    *nchoosek(N2,p)*(C0Coeff^N1)*(C1Coeff^(N2-p))*(C3Coeff^p)...
                    *sqrt(factorial(N1))*sqrt(factorial(N2-p))*sqrt(factorial(p));
                Ket42(k,1) = N1;
                Ket42(k,2) = N2;
                %KetN1_4(k,C0Pos+3) = N1;
                Ket42(k,C1Pos+4) = N2-p;
                Ket42(k,C3Pos+4) = p;
            end
            
        elseif N1 == 5;
            k=0;
            for p = 0:N2
                N1;
                N2;
                k=k+1;
                Ket51(k,3) = 1/(sqrt(factorial(N1))*sqrt(factorial(N2)))...
                    *nchoosek(N2,p)*(C0Coeff^N1)*(C1Coeff^(N2-p))*(C3Coeff^p)...
                    *sqrt(factorial(N1))*sqrt(factorial(N2-p))*sqrt(factorial(p));
                Ket51(k,1) = N1;
                Ket51(k,2) = N2;
                %KetN1_5(k,C0Pos+3) = N1;
                Ket51(k,C1Pos+4) = N2-p;
                Ket51(k,C3Pos+4) = p;
            end
            
        elseif N1 == 6;
            k=0;
            for p = 0:N2
                N1;
                N2;
                k=k+1;
                C0Coeff
                C0Pos
                C1Coeff
                C1Pos
                %C3Coeff
                %C3Pos
                Ket60(k,3) = 1/(sqrt(factorial(N1))*sqrt(factorial(N2)))...
                    *nchoosek(N2,p)*(C0Coeff^N1)*(C1Coeff^(N2-p))*(C3Coeff^p)...
                    *sqrt(factorial(N1))*sqrt(factorial(N2-p))*sqrt(factorial(p));
                Ket60(k,1) = N1;
                Ket60(k,2) = N2;
                Ket60(k,C0Pos+4) = N1;
                Ket60(k,C1Pos+4) = N2-p;
                Ket60(k,C3Pos+4) = p;
            end
            
            
            
            
            
            
        end
    end
  %% End Single BI
elseif C0Exist==1 && C2Exist==1 && C4Exist ==0 && C1Exist==1 && C3Exist==0
    
    %%%bixsingle
elseif C0Exist==1 && C2Exist==1 && C4Exist ==0 && C1Exist==1 && C3Exist==1
    %%%bixbi
elseif C0Exist==1 && C2Exist==1 && C4Exist ==1 && C1Exist==1 && C3Exist==0
    %%%trixsingle
elseif C0Exist==1 && C2Exist==1 && C4Exist ==1 && C1Exist==1 && C3Exist==1
    %%%trixbi
end

%%% finds cross of Tm with Wavefunction
aa=0;
[x,y]  =size(Ket06);
for CountA = 1:(x)
       for CountB= 1:length(Wavefunction)
        %KetN1_0(CountA,4:25)
        if Ket06(CountA,4:25) == Wavefunction(CountB,2:23)
            KetN = Ket06(CountA,1)
            Wv = Wavefunction(CountB,1)
            aa=aa+1;
            Ket06Val(aa,1) = (Ket06(CountA,3))*Wavefunction(CountB,1);
        end
    end
end

aa=0;
[x,y]  =size(Ket24);
for CountA = 1:(x)
       for CountB= 1:length(Wavefunction)
        
        if Ket24(CountA,4:25) == Wavefunction(CountB,2:23)
            KetN = Ket24(CountA,4:25);
            Wv = Wavefunction(CountB,2:23);
            aa=aa+1;
            Ket24Val(aa,1) = (Ket24(CountA,3))*Wavefunction(CountB,1);
        end
    end
end

aa=0;
[x,y]  =size(Ket42);
Ket42Val = zeros(1);
for CountA = 1:(x)
       for CountB= 1:length(Wavefunction)
        %KetN1_0(CountA,4:25)
        if Ket42(CountA,4:25) == Wavefunction(CountB,2:23)
            KetN = Ket42(CountA,4:25);
            Wv = Wavefunction(CountB,2:23);
            aa=aa+1;
            Ket42Val(aa,1) = (Ket42(CountA,3))*Wavefunction(CountB,1);
        end
    end
end

aa=0;
[x,y]  =size(Ket60);
for CountA = 1:(x)
       for CountB = 1:length(Wavefunction)
           
        %KetN1_0(CountA,4:25)
        if Ket60(CountA,4:25) == Wavefunction(CountB,2:23)
            KetN = Ket60(CountA,4:25)
            Wv = Wavefunction(CountB,1);
            aa=aa+1;
            Ket60Val(aa,1) = (Ket60(CountA,3))*Wavefunction(CountB,1);
            CountBPrint = CountB
        end
    end
end

Ket06Tot = (sum(Ket06Val.^2));
Ket24Tot = (sum(Ket24Val.^2));
Ket42Tot = (sum(Ket42Val.^2));
Ket60Tot = (sum(Ket60Val.^2));

Overlaps = [Ket60Tot Ket24Tot Ket42Tot Ket06Tot]
bar(Overlaps)
Omega
toc