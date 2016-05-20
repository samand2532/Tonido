clc; clear all; close all;
delete('UMatFinal.csv');delete('LMat.csv');delete('VMatFinalinA.csv');
AllBras = csvread('AllBras.csv');
AllKets = AllBras';
N = 6;
Lmax = 8;
A = 0.00;
g = 1;
MomValue = [-1 0 1 2 3 4 5 6 7 8 9 10 0]';

LmatFinal = zeros(length(AllBras));
UMatFinal = zeros(length(AllBras));
VMatFinal = zeros(length(AllBras));
%%%%%%N Term%%%%%%
NMatFinal = eye(length(AllBras)).*N;
%%%%%N Term End%%%%%%%

%%%%%% L Term%%%%%%%%
for BraL = 1:length(AllBras)
    for KetL = 1:length(AllKets)
        
        InitialBraL = AllBras(BraL,:);
        InitialKetL = AllKets(:,KetL);
        
        [OrthL] = Orthogonal(InitialBraL',InitialKetL);
        if OrthL == 1
            LMatFinal(BraL,KetL) = sum(InitialKetL.*MomValue);
        end
    end
end
csvwrite('LMat.csv',LMatFinal)
%%%%%%L Term End%%%%%%%%%%%%

%%%%%%%U Term%%%%%%%%
% generates m1:m4 possibilities
nn=0;qq=0;
for m1 = -1:Lmax
    for m2 = -1:Lmax
        for m3 = -1:Lmax
            for m4 = -1:Lmax
                nn=nn+1;
                StorageMatrixAll(nn,1) = m1;
                StorageMatrixAll(nn,2) = m2;
                StorageMatrixAll(nn,3) = m3;
                StorageMatrixAll(nn,4) = m4;
            end
        end
    end
end

for pp=1:(length(StorageMatrixAll))
    mL = StorageMatrixAll(pp,1)+StorageMatrixAll(pp,2);
    mR = StorageMatrixAll(pp,3)+StorageMatrixAll(pp,4);
    if (mL == mR) & (mL <= Lmax) & (mL ~=0)
        qq=qq+1;
        SortingMat(qq,:) = StorageMatrixAll(pp,:);
    end
end
AllPossibleArrangements = [0 0 0 0; SortingMat];
%section following removes all occurances of more than 1 lot of -1
remove = AllPossibleArrangements == -1;
remove = sum(remove,2);
jj=0;
for ii = 1:size(remove)
    if (remove(ii)>1)
    else jj = jj + 1; 
        AllPossibleArrangementsClean(jj,:) = AllPossibleArrangements(ii,:);
    end
end
APAC = AllPossibleArrangementsClean;

for KetsU = 1:length(AllKets)
    for BrasU = 1:length(AllBras)
        
        InitialBraU = AllBras(BrasU,:)';
        InitialKetU = AllKets(:,KetsU);
        
        %%% ann/crea vectors
        for zz = 1:length(APAC);
            
        Constm4Tmp=0;
        Constm3Tmp=0;
        Constm2Tmp=0;
        Constm1Tmp=0;
        m1 = APAC(zz,1)+2;
        m2 = APAC(zz,2)+2;
        m3 = APAC(zz,3)+2;
        m4 = APAC(zz,4)+2;
        m1Vect = zeros(13,1);
        m2Vect = zeros(13,1);
        m3Vect = zeros(13,1);
        m4Vect = zeros(13,1);
        m1Vect(m1,1) = 1;
        m2Vect(m2,1) = 1;
        m3Vect(m3,1) = -1;
        m4Vect(m4,1) = -1;
        
        m4Const = sqrt(InitialKetU(m4,1));
        m4trans = InitialKetU + m4Vect;
        m3Const = sqrt(m4trans(m3,1));
        m3trans = m4trans + m3Vect;
        m2Const = sqrt(m3trans(m2,1)+1);
        m2trans = m3trans + m2Vect;
        m1Const = sqrt(m2trans(m1,1)+1);
        FinalKetU = m2trans + m1Vect;
        OperatorConst = m4Const*m3Const*m2Const*m1Const;
        [OrthU] = Orthogonal(InitialBraU,FinalKetU);
       
        if OrthU == 1
            I2 = factorial( (abs(m1-2)+abs(m2-2)+abs(m3-2)+abs(m4-2))/2 );
            PiTerm = sqrt((factorial(InitialKetU(13,1))*(factorial(InitialKetU(13,1)))*(factorial(InitialKetU(13,1)))*(factorial(InitialKetU(13,1)))) / (factorial(InitialKetU(13,1) + abs(m1-2))*factorial(InitialKetU(13,1) + abs(m2-2))*factorial(InitialKetU(13,1) + abs(m3-2))*factorial(InitialKetU(13,1) + abs(m4-2))));
            MomTerm = 
            Total=I2*PiTerm*OperatorConst
        end
        end
    end
end

        
        

