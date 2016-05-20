clc; clear all;

Nt = 6;

%%%% CHANGE THESE VALUES ONLY!!!!%%%%
x = 0.98948;
x1 = -0.14224;
Eig = struct2array(load('V0_77725.mat'));
%Eig = csvread('V_077721759.csv');
%%%%%%%%%%%%%%%%%%%%%%

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
%%%%%momentum basis used initially, next part finds fidelity with original
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
