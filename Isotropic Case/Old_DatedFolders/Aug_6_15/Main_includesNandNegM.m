clc; clear all; close all;
delete('UMatFinal.csv'); delete('LMat.csv'); delete('VMatFinalinA.csv');
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
        
        %[OrthL] = Orthogonal(InitialBraL',InitialKetL);
        %if OrthL == 1
        if InitialBraL' == InitialKetL
            LmatFinal(BraL,KetL) = sum(InitialKetL.*MomValue);
        end
    end
end
csvwrite('LMat.csv',LmatFinal)
%%%%%%L Term End%%%%%%%%%%%%

%%%%%%%U Term%%%%%%%%
% generates m1:m4 possibilities

nn=0;qq=0;
StorageMatrixAll = zeros((Lmax+2)^4*5,5);
for m1 = -1:Lmax
    for m2 = -1:Lmax
        for m3 = -1:Lmax
            for m4 = -1:Lmax
                for n_i = 0:4
                nn=nn+1;
                StorageMatrixAll(nn,1) = m1;
                StorageMatrixAll(nn,2) = m2;
                StorageMatrixAll(nn,3) = m3;
                StorageMatrixAll(nn,4) = m4;
                StorageMatrixAll(nn,5) = n_i;
            end
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
AllPossibleArrangements = [0 0 0 0 0;0 0 0 0 1;0 0 0 0 2;0 0 0 0 3;0 0 0 0 4; SortingMat];
%section following removes all occurances of more than 1 lot of -1
remove = AllPossibleArrangements == -1;
remove = sum(remove,2);
jj=0;
for ii = 1:size(remove)
    if (remove(ii)>1)
    else jj = jj + 1; 
        AllPossibleArrangements1(jj,:) = AllPossibleArrangements(ii,:);
    end
end
jj1 = 0;
for jj=1:length(AllPossibleArrangements1)
    
    if ((AllPossibleArrangements1(jj,1) ==-1) || (AllPossibleArrangements1(jj,2) ==-1) || (AllPossibleArrangements1(jj,3) ==-1) || (AllPossibleArrangements1(jj,4) ==-1)) && ((AllPossibleArrangements1(jj,5) > 0))
    
    else jj1 = jj1+1; 
        AllPossibleArrangementsClean(jj1,:) = AllPossibleArrangements1(jj,:);
    end
end
APAC = AllPossibleArrangementsClean;

        m1Vect = zeros(13,1);
        m2Vect = zeros(13,1);
        m3Vect = zeros(13,1);
        m4Vect = zeros(13,1);

for KetsU = 1:length(AllKets)
    for BrasU = 1:length(AllBras)
        
        InitialBraU = AllBras(BrasU,:)';
        InitialKetU = AllKets(:,KetsU);
        
        %%% ann/crea vectors
        for zz = 1:length(APAC);
            
        m1 = APAC(zz,1)+2;
        m2 = APAC(zz,2)+2;
        m3 = APAC(zz,3)+2;
        m4 = APAC(zz,4)+2;
        m1Vect = m1Vect*0;
        m2Vect = m2Vect*0;
        m3Vect = m3Vect*0;
        m4Vect = m4Vect*0;
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
        
       % [OrthU] = Orthogonal(InitialBraU,FinalKetU);
       
       % if OrthU == 1
       if InitialBraU == FinalKetU
           
           OperatorConst = m4Const*m3Const*m2Const*m1Const;
           I2 = factorial( (abs(m1-2)+abs(m2-2)+abs(m3-2)+abs(m4-2))/2 );
           
           if APAC(zz,5) == 0
               nk1 = 0; nk2 = 0; nk3 = 0; nk4 = 0;
           elseif APAC(zz,5) == 1
               nk1 = 1; nk2 = 0; nk3 = 0; nk4 = 0;
           elseif APAC(zz,5) == 2
               nk1 = 0; nk2 = 1; nk3 = 0; nk4 = 0;
           elseif APAC(zz,5) == 3
               nk1 = 0; nk2 = 0; nk3 = 1; nk4 = 0;
           elseif APAC(zz,5) == 4
               nk1 = 0; nk2 = 0; nk3 = 0; nk4 = 1;
           end
               
               PiTerm = sqrt((1) / (factorial(nk1 + abs(m1-2))*factorial(nk2 + abs(m2-2))*factorial(nk3 + abs(m3-2))*factorial(nk4 + abs(m4-2))));
               
               MomTerm = 1/(2^(0.5*LmatFinal(BrasU,KetsU)));
               TotalU=(g/(4*pi))*(MomTerm)*(I2)*(PiTerm)*(OperatorConst);
               UMatFinal(BrasU,KetsU) = UMatFinal(BrasU,KetsU) + TotalU;
           end
       end
    end
end
csvwrite('UMatFinal.csv',UMatFinal);
tt=0;
for Omega = 0.6:0.01:0.85;
    tt=tt+1;
    NAll = NMatFinal;
    UAll = csvread('UMatFinal.csv');
    LAllread = csvread('LMat.csv');
    LAll = LAllread.*(1-Omega);
    Total = NAll + UAll + LAll;
    Eig = eig(Total);
    
    OmegaIncreaseMatrix(1,tt) = Omega;
    OmegaEnergyMatrix(:,tt) = Eig;
    OmegaMatRestricted = OmegaEnergyMatrix(1:9,:);
    
end
        
plot(OmegaIncreaseMatrix, OmegaEnergyMatrix)
xlabel('Omega')
ylabel('<E>')
%ylim([8 9.2])
%xlim([0.6 0.85])
hold off;
        

