clc; clear all;
%V = struct2array(load('EigenMat0796.mat'));
%load('V0.8.mat')
%load('EigenMat0796.mat')
AllBras = csvread('StatesLLL.csv');
AllKets = AllBras';
%V = struct2array(load('V_077721759.csv'));
V = csvread('V_077721759.csv');
eigenvectors = V(:,1);
format long
eigMat = zeros(9);

creaVect = zeros(19,1);
annVect = zeros(19,1);

for Kets = 1:39;
for k = [2 4 6 8 10 12 14 16 18];
for l = [2 4 6 8 10 12 14 16 18];
    C1 = eigenvectors(Kets,1);
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
        C2 = eigenvectors(a,1);
        if Bra' == FinalKet
            eigMat(k/2,l/2) = eigMat(k/2,l/2) + (creaConst*annConst*C1*C2);
        end
    end
end
end
end

eigMat;
[V1, D] = eig(eigMat);
dia = diag(eigMat);
N = 1/sqrt(sum(dia).^2);
Norm = dia.*N;
SumNormEigVal = sum(Norm );

D;

V1;