%%% Combined Mk II
clear all; clc;
delete('UMatFinal.csv');delete('LMat.csv');delete('VMatFinalincA.csv')
AllBras = csvread('AllBras.csv');
AllKets = AllBras';
N=6;
Lmax=8;
A=0.0;
g=1;
MomValue = [-1 0 1 2 3 4 5 6 7 8 9 10 0]'; % first 2 0's are from n and m, dont cont. to LmatMatrix
%%%%%%%%%%%%%%NTerm%%%%%%%%%%%%%%%%%%
NMatFinal = eye(length(AllBras)).*N;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%LTERM&&&&&&&&&&&&&&&&&
LMatFinal = zeros(length(AllBras));
UMatFinal = zeros(length(AllBras));
VMatFinal = zeros(length(AllBras));
for BraL = 1:length(AllBras)
    for KetL = 1:length(AllKets)
        
        InitialBraL = AllBras(BraL,:);
        InitialKetL = AllKets(:,KetL);
        
        [OrthL] = Orthogonal(InitialBraL', InitialKetL);
        if OrthL == 1
            LMatFinal(BraL,KetL) = sum(InitialKetL.*MomValue);
        end
    end
end
csvwrite('LMat.csv',LMatFinal)

%% UTerm
%%%%%%%%%%%%%START U TERM%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%% GENERATES ALL POSSIBLE m arrangements upto Mt = 8.
%%%%% Generates all arrangements of m1:m4.
n=0;q=0;
for m1 = 0:(Lmax+1);
    for m2 = 0:(Lmax+1);
        for m3 = 0:(Lmax+1);
            for m4 = 0:(Lmax+1);
                n = n+1;
                StorageMatrixAll(n,1)=m1-1;
                StorageMatrixAll(n,2)=m2-1;
                StorageMatrixAll(n,3)=m3-1;
                StorageMatrixAll(n,4)=m4-1;
            end
        end
    end
end

for p=1:(length(StorageMatrixAll))
    mL = StorageMatrixAll(p,1)+StorageMatrixAll(p,2);
    mR = StorageMatrixAll(p,3)+StorageMatrixAll(p,4);
    if (mL == mR & mL <= (Lmax) & mL ~= 0)
        q=q+1;
        c(q,:)=StorageMatrixAll(p,:);
    end 
end
%%%%% This section inserts 0 0 0 0 once.
AllPossibleArrangements = [0 0 0 0; c];

remove = AllPossibleArrangements==-1; % remove is the boolean matrix of APA that is equal to -1
remove = sum(remove,2); % this sums it to over the row wto see if there  are more than 1 lots of (-1)
j=0;
for i=1:size(remove)
    if (remove(i)>1)
    else j=j+1; APA_2(j,:)=AllPossibleArrangements(i,:);
    end
end

%%%%%%%%%%%% Big constant time%%%%%%%
%BigConst = (factorial(m1 + m2)/ sqrt(factorial(m1)*factorial(m2)*factorial(m3)*factorial(m4))) * (1/ 2^(m1+m2+1))*(g/p);
for q=1:length(AllPossibleArrangements)
    AllPossibleArrangements(q,5) = (factorial(abs(AllPossibleArrangements(q,1)) + abs(AllPossibleArrangements(q,2)))/ sqrt(factorial(abs(AllPossibleArrangements(q,1)))*factorial(abs(AllPossibleArrangements(q,2)))*factorial(abs(AllPossibleArrangements(q,3)))*factorial(abs(AllPossibleArrangements(q,4))))) * (1/ 2^(abs(AllPossibleArrangements(q,1))+abs(AllPossibleArrangements(q,2))+1))*(g/pi);
end
%%%%%%AllPossibleArrangements Has Form |m1|m2|m3|m4|BigConstant
%                                      |  |  |  |  |
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for KetsU=1:length(AllKets)
    for BrasU=1:length(AllBras)
        
        InitialBraU = AllBras(BrasU,:)';
        InitialKetU = AllKets(:,KetsU);
        
        %%%%%%%%%%%%%%% m ann / crea vectors
        for z=1:length(AllPossibleArrangements);
            Constm4Tmp=0;
            Constm3Tmp=0;
            Constm4=0;
            Constm3=0;
            Constm2=0;
            Constm1=0;
            
            m1 = AllPossibleArrangements(z,1);
            m2 = AllPossibleArrangements(z,2); 
            m3 = AllPossibleArrangements(z,3); 
            m4 = AllPossibleArrangements(z,4);
            m1Vect = zeros(13,1); 
            m2Vect = zeros(13,1); 
            m3Vect = zeros(13,1); 
            m4Vect = zeros(13,1);
            m1Vect((m1+2),1)=1;
            m2Vect((m2+2),1)=1;
            m3Vect((m3+2),1)=-1;
            m4Vect((m4+2),1)=-1;
            
            m4Const = sqrt(InitialKetU((m4+2),1));
            m4trans = InitialKetU + m4Vect;
            m3Const = sqrt(m4trans((m3+2),1));
            m3trans = m4trans + m3Vect;
            m2Const = sqrt(m3trans((m2+2),1)+1);
            m2trans = m3trans + m2Vect;
            m1Const = sqrt(m2trans((m1+2),1)+1);
            FinalKetU = m2trans + m1Vect;
            
            %%%%%%%%Orthoganality Test
            [OrthU] = Orthogonal(InitialBraU,FinalKetU);
            %%%%%%%%If orthogonal find constants and save to matrix
            
            ConstComb = OrthU*m4Const*m3Const*m2Const*m1Const;
            Matelement1(z,1) = ConstComb*AllPossibleArrangements(z,5);
            
        end
        UMatFinal(BrasU,KetsU) = 0.5*sum(Matelement1);
    end
end
csvwrite('UMatFinal.csv',UMatFinal)

%%%%%%%%VTERM HERE
AllKets =  AllBras';
VFinal = zeros(124);


for bra = 1:length(AllBras)
    for ket = 1:length(AllKets)
        
        InitialKet = AllKets(:,ket);
        InitialBra = AllBras(bra,:)';
        
        %m0
        P2m0Const = sqrt(2) * sqrt(InitialKet(1,1)+1) * sqrt(InitialKet(3,1));
        P2m0Ket = InitialKet + [  1 0 0 0 0 0 0 0 0 0 0 0 0]' + [  0 0 -1 0 0 0 0 0 0 0 0 0 0]';
        
        %m1
        P2m1Const = sqrt(6) * sqrt(InitialKet(2,1)+1) * sqrt(InitialKet(4,1));
        P2m1Ket = InitialKet + [  0 1 0 0 0 0 0 0 0 0 0 0 0]' + [  0 0 0 -1 0 0 0 0 0 0 0 0 0]';
        
        %m2
        P1m2Const = sqrt(2) * sqrt(InitialKet(3,1)+1) * sqrt(InitialKet(1,1));
        P1m2Ket = InitialKet + [  0 0 1 0 0 0 0 0 0 0 0 0 0]' + [  -1 0 0 0 0 0 0 0 0 0 0 0 0]';
        P2m2Const = sqrt(12) * sqrt(InitialKet(3,1)+1) * sqrt(InitialKet(5,1));
        P2m2Ket = InitialKet + [  0 0 1 0 0 0 0 0 0 0 0 0 0]' + [  0 0 0 0 -1 0 0 0 0 0 0 0 0]';
        %m3
        P1m3Const = sqrt(6) * sqrt(InitialKet(4,1)+1) * sqrt(InitialKet(2,1));
        P1m3Ket = InitialKet + [  0 0 0 1 0 0 0 0 0 0 0 0 0]' + [  0 -1 0 0 0 0 0 0 0 0 0 0 0]';
        P2m3Const = sqrt(20) * sqrt(InitialKet(4,1)+1) * sqrt(InitialKet(6,1));
        P2m3Ket = InitialKet + [  0 0 0 1 0 0 0 0 0 0 0 0 0]' + [  0 0 0 0 0 -1 0 0 0 0 0 0 0]';
        %m4
        P1m4Const = sqrt(12) * sqrt(InitialKet(5,1)+1) * sqrt(InitialKet(3,1));
        P1m4Ket = InitialKet + [  0 0 0 0 1 0 0 0 0 0 0 0 0]' + [  0 0 -1 0 0 0 0 0 0 0 0 0 0]';
        P2m4Const = sqrt(30) * sqrt(InitialKet(5,1)+1) * sqrt(InitialKet(7,1));
        P2m4Ket = InitialKet + [  0 0 0 0 1 0 0 0 0 0 0 0 0]' + [  0 0 0 0 0 0 -1 0 0 0 0 0 0]';
        %m5
        P1m5Const = sqrt(20) * sqrt(InitialKet(6,1)+1) * sqrt(InitialKet(4,1));
        P1m5Ket = InitialKet + [  0 0 0 0 0 1 0 0 0 0 0 0 0]' + [  0 0 0 -1 0 0 0 0 0 0 0 0 0]';
        P2m5Const = sqrt(42) * sqrt(InitialKet(6,1)+1) * sqrt(InitialKet(8,1));
        P2m5Ket = InitialKet + [  0 0 0 0 0 1 0 0 0 0 0 0 0]' + [  0 0 0 0 0 0 0 -1 0 0 0 0 0]';
        %m6
        P1m6Const = sqrt(30) * sqrt(InitialKet(7,1)+1) * sqrt(InitialKet(5,1));
        P1m6Ket = InitialKet + [  0 0 0 0 0 0 1 0 0 0 0 0 0]' + [  0 0 0 0 -1 0 0 0 0 0 0 0 0]';
        P2m6Const = sqrt(56) * sqrt(InitialKet(7,1)+1) * sqrt(InitialKet(9,1));
        P2m6Ket = InitialKet + [  0 0 0 0 0 0 1 0 0 0 0 0 0]' + [  0 0 0 0 0 0 0 0 -1 0 0 0 0]';
        %m7
        P1m7Const = sqrt(42) * sqrt(InitialKet(8,1)+1) * sqrt(InitialKet(6,1));
        P1m7Ket = InitialKet + [  0 0 0 0 0 0 0 1 0 0 0 0 0]' + [  0 0 0 0 0 -1 0 0 0 0 0 0 0]';
        P2m7Const = sqrt(72) * sqrt(InitialKet(8,1)+1) * sqrt(InitialKet(10,1));
        P2m7Ket = InitialKet + [  0 0 0 0 0 0 0 1 0 0 0 0 0]' + [  0 0 0 0 0 0 0 0 0 -1 0 0 0]';
        %m8
        P1m8Const = sqrt(56) * sqrt(InitialKet(9,1)+1) * sqrt(InitialKet(7,1));
        P1m8Ket = InitialKet + [  0 0 0 0 0 0 0 0 1 0 0 0 0]' + [  0 0 0 0 0 0 -1 0 0 0 0 0 0]';
        P2m8Const = sqrt(90) * sqrt(InitialKet(9,1)+1) * sqrt(InitialKet(11,1));
        P2m8Ket = InitialKet + [  0 0 0 0 0 0 0 0 1 0 0 0 0]' + [  0 0 0 0 0 0 0 0 0 0 -1 0 0]';
        
        
        
        [P2m0Ortho] = Orthogonal (P2m0Ket, InitialBra);
        [P2m1Ortho] = Orthogonal (P2m1Ket, InitialBra);
        [P1m2Ortho] = Orthogonal (P1m2Ket, InitialBra);
        [P2m2Ortho] = Orthogonal (P2m2Ket, InitialBra);
        [P1m3Ortho] = Orthogonal (P1m3Ket, InitialBra);
        [P2m3Ortho] = Orthogonal (P2m3Ket, InitialBra);
        [P1m4Ortho] = Orthogonal (P1m4Ket, InitialBra);
        [P2m4Ortho] = Orthogonal (P2m4Ket, InitialBra);
        [P1m5Ortho] = Orthogonal (P1m5Ket, InitialBra);
        [P2m5Ortho] = Orthogonal (P2m5Ket, InitialBra);
        [P1m6Ortho] = Orthogonal (P1m6Ket, InitialBra);
        [P2m6Ortho] = Orthogonal (P2m6Ket, InitialBra);
        [P1m7Ortho] = Orthogonal (P1m7Ket, InitialBra);
        [P2m7Ortho] = Orthogonal (P2m7Ket, InitialBra);
        [P1m8Ortho] = Orthogonal (P1m8Ket, InitialBra);
        [P2m8Ortho] = Orthogonal (P2m8Ket, InitialBra);
        
        if P2m0Ortho ==1
            VFinal(bra,ket) = VFinal(bra,ket) + P2m0Const;
        end
        if P2m1Ortho ==1
            VFinal(bra,ket) = VFinal(bra,ket) + P2m1Const;
        end
        if P1m2Ortho ==1
            VFinal(bra,ket) = VFinal(bra,ket) + P1m2Const;
        end
        if P2m2Ortho ==1
            VFinal(bra,ket) = VFinal(bra,ket) + P2m2Const;
        end
        if P1m3Ortho ==1
            VFinal(bra,ket) = VFinal(bra,ket) + P1m3Const;
        end
        if P2m3Ortho ==1
            VFinal(bra,ket) = VFinal(bra,ket) + P2m3Const;
        end
        if P1m4Ortho ==1
            VFinal(bra,ket) = VFinal(bra,ket) + P1m4Const;
        end
        if P2m4Ortho ==1
            VFinal(bra,ket) = VFinal(bra,ket) + P2m4Const;
        end
        if P1m5Ortho ==1
            VFinal(bra,ket) = VFinal(bra,ket) + P1m5Const;
        end
        if P2m5Ortho ==1
            VFinal(bra,ket) = VFinal(bra,ket) + P2m5Const;
        end
        if P1m6Ortho ==1
            VFinal(bra,ket) = VFinal(bra,ket) + P1m6Const;
        end
        if P2m6Ortho ==1
            VFinal(bra,ket) = VFinal(bra,ket) + P2m6Const;
        end
        if P1m7Ortho ==1
            VFinal(bra,ket) = VFinal(bra,ket) + P1m7Const;
        end
        if P2m7Ortho ==1
            VFinal(bra,ket) = VFinal(bra,ket) + P2m7Const;
        end
        if P1m8Ortho ==1
            VFinal(bra,ket) = VFinal(bra,ket) + P1m8Const;
        end
        if P2m8Ortho ==1
            VFinal(bra,ket) = VFinal(bra,ket) + P2m8Const;
        end
    end
end

VFinalincA = VFinal.*A;
csvwrite('VMatFinalincA.csv',VFinalincA)

%%%%%%%%%%%%%%%%%


%% Combining
%For varying Omega
for Omega = 0.6:0.01:0.85;
    n=n+1;
    
    NAll = NMatFinal;
    UAll = csvread('UMatFinal.csv');
    LAllread = csvread('LMat.csv');
    LAll = LAllread.*(1-Omega);
    VAll = csvread('VMatFinalincA.csv');
    Total = NAll +UAll+(VAll)+(LAll);
    Eig = eig(Total);
    
    
    OmegaIncreaseMatrix(1,n) = Omega;
    OmegaEnergyMatrix(:,n) = Eig;
    %%% Following line selects only the lowest 8 eigenstates
    OmegaMatRestricted = OmegaEnergyMatrix(1:9,:);
     
end

plot(OmegaIncreaseMatrix,OmegaEnergyMatrix);
xlabel('Omega')
ylabel('<E>')
ylim([8 9.2])
xlim([0.6 0.85])
hold off;

