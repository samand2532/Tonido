clc; clear all;
eigenstates = load('V0_777216.mat');
%load('V0_777216.mat')
%eigenstatesAll = struct2array(load('TESTEigenMatrix1.mat'));
%load('TESTEigenMatrixOmega.mat')

AllBras = csvread('StatesLLL.csv');
AllKets = AllBras';
tic

eigMat = zeros(9);

creaVect = zeros(19,1);
annVect = zeros(19,1);
count = 0;
for qw = 1:481;
    count = count+1;
    if count == 1000
        count
    end
    if count == 2000
        count
    end
    
    if count == 3000 
        count
    end
    if count == 4000
        count
    end
      if count == 5000
        count
    end
  %  eigenstates = eigenstatesAll(:,qw);
    eigMat = eigMat.*0;

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
    dia = diag(eigMat);
sumdiag = sum(dia);
diaStorMat(:,count) = dia;

end
end
end
end
toc
plot(TESTEigenMatrixOmega,diaStorMat(1:3,:))
legend('Psi_0','Psi_1','Psi_2')
