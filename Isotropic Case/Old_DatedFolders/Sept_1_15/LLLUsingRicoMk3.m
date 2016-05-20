clear all; clc;


if exist(('LMat.csv'),'file') == 2
    delete ('LMat.csv')
end
if exist(('UMatFinal.csv'),'file') == 2
    delete ('LMat.csv')
end

    
AllBras = csvread('AllBrasLLL.csv');
AllKets = AllBras';
Lmax = 8;
g = 1;
A = 0.03;
N = 6;
MomValue = [ -1 0 1 2 3 4 5 6 7 8 9 10 0]';
%NoOfVectsLookedAt = 3; %FOR DEBUGGING!  SAVES HAVING TO CALC FOR ALL VECTORS< AND CHANGING INDIV VALUES
Nmat = zeros(39);
Lmat = zeros(39);
UmatNoConst = zeros(39);
Vmat = zeros(39);

%%%%%%NTerm Follows
NMatFinal = eye(39).*N;
%%%%%%%NTermEnd%%%%

%%%%%L Term Begin%%%%
for BraL = 1:length(AllBras)
    for KetL = 1:length(AllBras)
        InitialBraL = AllBras(BraL,:);
        InitialKetL = AllKets(:,KetL);       
        if InitialBraL == InitialKetL'
            LMatFinal(BraL,KetL) = dot(InitialKetL,MomValue);
        end
    end
end
csvwrite('LMat.csv', LMatFinal);
%%%%L Term End%%%
%%%%UtermStart
nn=0;qq=0;rr=0;
StorageMatrixAll = zeros(4^Lmax,5);

for m1 = 1:(Lmax+1);
    for m2 = 1:(Lmax+1);
        for m3 = 1:(Lmax+1);
            for m4 = 1:(Lmax+1);
                nn = nn+1;
                StorageMatrixAll(nn,1) = m1-1;
                StorageMatrixAll(nn,2) = m2-1;
                StorageMatrixAll(nn,3) = m3-1;
                StorageMatrixAll(nn,4) = m4-1;
            end
        end
    end
end

for pp = 1:(4^Lmax)
    mL = StorageMatrixAll(pp,1)+StorageMatrixAll(pp,2);
    mR = StorageMatrixAll(pp,3)+StorageMatrixAll(pp,4);
    if (mL == mR & mL > 0 & mL <= Lmax)
        qq = qq + 1;
        SM(qq,:) = StorageMatrixAll(pp,:);
    end
end
APA = [0 0 0 0 0; SM];

m1Vect = zeros(13,1);
m2Vect = zeros(13,1);
m3Vect = zeros(13,1);
m4Vect = zeros(13,1);
ee=0;
MomTest=zeros(200,1);
for BrasU = 1:length(AllKets)
    for KetsU = 1:length(AllBras)
        InitialBraU = AllBras(BrasU,:)';
        InitialKetU = AllKets(:,KetsU);
        
        for rr = 1:length(APA);
            %%%% NOT SURE IF +2 or what here :::((((((
            m1 = APA(rr,1)+2;
            m2 = APA(rr,2)+2;
            m3 = APA(rr,3)+2;
            m4 = APA(rr,4)+2;
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
            
            if FinalKetU == InitialBraU
                OpConst = m1Const*m2Const*m3Const*m4Const
                I2 = gamma((abs(m1)-2 + abs(m2)-2 + abs(m3)-2 + abs(m4)-2)/2 + 1)
                
                if APA(rr,5) == 0
                    nk1 = 0; nk2 = 0; nk3 = 0; nk4 = 0;
                elseif APA(rr,5) == 1
                    nk1 = 1; nk2 = 0; nk3 = 0; nk4 = 0;
                elseif APA(rr,5) == 2
                    nk1 = 0; nk2 = 1; nk3 = 0; nk4 = 0;
                elseif APA(rr,5) == 3
                    nk1 = 0; nk2 = 0; nk3 = 1; nk4 = 0;
                elseif APA(rr,5) == 4
                    nk1 = 0; nk2 = 0; nk3 = 0; nk4 = 1;
                end
                
                PiTerm = (factorial(nk1)*factorial(nk2)*factorial(nk3)*factorial(nk4))/...
                    (factorial(nk1+abs(m1)-2)*factorial(nk2+abs(m2)-2)*factorial(nk3+abs(m3)-2)*factorial(nk4+abs(m4)-2));
                MomTerm = 1/(2^((dot(InitialBraU,abs(MomValue)))/2))
                UMid = OpConst*I2*PiTerm*MomTerm;
                UmatNoConst(BrasU,KetsU) = UmatNoConst(BrasU,KetsU) + UMid;
            end
        end
    end
end
UMatFinal = (g/4*pi).* UmatNoConst;
csvwrite('UMatFinal.csv', UMatFinal);       
%%%End of U Term
%%%%Start of V term



%%%%% ENd of V term

%% VARYING OMEGA
ss=0;
for Omega = 0.6:0.001:1.5
    ss=ss+1;
    
    NAll = NMatFinal;
    UAll = csvread('UMatFinal.csv');
    LAll = (csvread('LMat.csv')).*(1-Omega);
    Total = NAll+UAll+LAll;
    Eig = eig(Total);
    
    OmegaIncreaseMatrix(1,ss) = Omega;
    OmegaEnergyMatrix(:,ss) = Eig;
    %OmegaMatRestricted = OmegaEnergyMatrix(1:9,:);
    
end

plot(OmegaIncreaseMatrix,OmegaEnergyMatrix);
xlabel('Omega')
ylabel('<E>')
ylim([8 9.2])
xlim([0.6 0.85])
hold off;
    



                
                    
                
        
                
                
                
                
                
                
                
                
                
                
                
                
                
                
