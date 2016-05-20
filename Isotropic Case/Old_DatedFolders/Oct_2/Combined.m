clc; clear all;

KPos = [ 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19; %Matlab position
         0 0 1 0 1 0 1 0 1 0  1  0  1  0  1  0  1  0  1 ; %n value
        -1 0 0 1 1 2 2 3 3 4  4  5  5  6  6  7  7  8  8]; %Mt value

MomValue = [ -1 0 0 1 1 2 2 3 3 4 4 5 5 6 6 7 7 8 8];

g = 1; N = 6; A = 0.015; LLVal = 1; % 1 = LLL, 2 = LL2
format short
tic
if LLVal == 1
    
    AllBras = csvread('StatesLLL.csv');
    AllKets = AllBras';
    
       
    NMat = eye(length(AllBras)).*N;
    for BraL = 1:length(AllBras)
        for KetL = 1:length(AllKets)
            InitialBraL = AllBras(BraL,:);
            InitialKetL = AllKets(:,KetL);
            
            if InitialBraL == InitialKetL'
                LMat(BraL,KetL) = dot(InitialKetL, MomValue);
            end
            
        end
    end
    
    deltaINT = csvread('LagIntLLL.csv');
    deltaPOT = csvread('LagPotLLL.csv');
    UMatNogConst = zeros(length(AllKets));
    
    l1Vect = zeros(size(InitialBraL))';
    l2Vect = zeros(size(InitialBraL))';
    k1Vect = zeros(size(InitialBraL))';
    k2Vect = zeros(size(InitialBraL))';
    %%%% INT PART
    for BrasU = 1:length(AllBras)
        BrasU
        for KetsU = 1:length(AllBras)
            InitialBraU = AllBras(BrasU,:);
            InitialKetU = AllKets(:,KetsU);
            
    for deltaCountINT = 1:length(deltaINT)
        k1 = deltaINT(deltaCountINT,1);
        k2 = deltaINT(deltaCountINT,2);
        l1 = deltaINT(deltaCountINT,3);
        l2 = deltaINT(deltaCountINT,4);
        nk1 = deltaINT(deltaCountINT,5);
        nk2 = deltaINT(deltaCountINT,6);
        n11 = deltaINT(deltaCountINT,7);
        nl2 = deltaINT(deltaCountINT,8);
        mk1 = deltaINT(deltaCountINT,9);
        mk2 = deltaINT(deltaCountINT,10);
        ml1 = deltaINT(deltaCountINT,11);
        ml2 = deltaINT(deltaCountINT,12);
        
        l1Vect = l1Vect.*0;
        l2Vect = l2Vect.*0;
        k1Vect = k1Vect.*0;
        k2Vect = k2Vect.*0;
        
        l1Vect(l1,1) = -1;
        l2Vect(l2,1) = -1;
        k1Vect(k1,1) = 1;
        k2Vect(k2,1) = 1;
        
        l2Const = sqrt(InitialKetU(l2,1));
        l2trans  = InitialKetU + l2Vect;
        l1Const = sqrt(l2trans(l1,1));
        l1trans = l2trans + l1Vect;
        k2Const = sqrt(l1trans(k2,1)+1);
        k2trans = l1trans + k2Vect;
        k1Const = sqrt(k2trans(k1,1)+1);
        FinalKet = k2trans + k1Vect;
      
        if FinalKet == InitialBraU'
            
            OpConst = l2Const*l1Const*k2Const*k1Const;
            OtherConst = deltaINT(deltaCountINT,13);  
            UMatNogConst(BrasU,KetsU) = UMatNogConst(BrasU,KetsU) + (OtherConst*OpConst);
        end
    end
        end
    end
    %%%%% POT TERM
    
    k1Vect = zeros(size(InitialBraL))';
    k2Vect = zeros(size(InitialBraL))';
    
    VMatNoA = zeros(length(AllKets));
    
    for BrasV = 1:length(AllBras)
        for KetsV = 1:length(AllBras)
            InitialBraV = AllBras(BrasV,:);
            InitialKetV = AllKets(:,KetsV);
    
            for deltaCountPOT = 1:length(deltaPOT)
                
                k1 = deltaPOT(deltaCountPOT,1);
                k2 = deltaPOT(deltaCountPOT,2);
                nk1 = deltaPOT(deltaCountPOT,3);
                nk2 = deltaPOT(deltaCountPOT,4);
                mk1 = deltaPOT(deltaCountPOT,5);
                mk2 = deltaPOT(deltaCountPOT,6);
                
                k1Vect = k1Vect.*0;
                k2Vect = k2Vect.*0;
                
                k1Vect(k1,1) = 1;
                k2Vect(k2,1) = -1;
                
                k2Const = sqrt(InitialKetV(k2,1));
                k2trans = InitialKetV + k2Vect;
                k1Const = sqrt(k2trans(k1,1)+1);
                FinalKetV = k2trans + k1Vect;
                
                if FinalKetV == InitialBraV'
                    
                    OpConst = k2Const*k1Const;
                    OtherConst = deltaPOT(deltaCountPOT,7);
                    %I1 = factorial((abs(mk1)+abs(mk2)+2)/2);
                    %PiT = sqrt(1/(factorial(abs(mk1))*factorial(abs(mk2))));
                    %OtherConst = I1*PiT;
                    VMatNoA(BrasV,KetsV) = VMatNoA(BrasV,KetsV) + (OtherConst*OpConst);
                end
            end
        end
    end
    
    
    
    n = 0;
for Omega = 0.777216;%0.6:0.0001:1.2;
    n=n+1;
    
    NAll = NMat;
    LAll = LMat*(1-Omega);
    UAll = UMatNogConst.*(g/(4*pi));
    UAlla=round(UAll,8);
    VMat = VMatNoA.*A;
    %VMat = round(VMata,8);
    Total = NAll +UAlla+(LAll)+VMat;
    Eig = eig(Total);
    [V D] = eig(Total);
    
    
    OmegaIncreaseMatrix(1,n) = Omega;
    OmegaEnergyMatrix(:,n) = Eig;
    %%% Following line selects only the lowest 8 eigenstates
    OmegaMatRestricted = OmegaEnergyMatrix(1:9,:);
    
    
end

plot(OmegaIncreaseMatrix,OmegaEnergyMatrix);
xlabel('Omega')
ylabel('<E>')
ylim([8.2 8.9])
xlim([0.7 0.8])
    
toc    
end

    
  
if LLVal == 2
    
    AllBras = csvread('StatesLLL_and_LL2.csv');
    AllKets = AllBras';
    
    NMat = eye(length(AllBras)).*N;
    for BraL = 1:length(AllBras)
        for KetL = 1:length(AllKets)
            InitialBraL = AllBras(BraL,:);
            InitialKetL = AllKets(:,KetL);
            
            if InitialBraL == InitialKetL'
                LMat(BraL,KetL) = dot(InitialKetL, MomValue);
            end
            
        end
    end
    
    %deltaINT = csvread('LagIntLL2.csv');
    %deltaINT = csvread('DeltaIntTry2LL2.csv');%%%% THIS IS TRY 2
    %deltaINT = csvread('DeltaIntTry3LL2.csv')
    deltaINT = csvread('DeltaIntLL2Try4.csv')
    deltaPOT = csvread('LagPotLL2.csv');
    UMatNogConst = zeros(length(AllKets));
    
    l1Vect = zeros(size(InitialBraL))';
    l2Vect = zeros(size(InitialBraL))';
    k1Vect = zeros(size(InitialBraL))';
    k2Vect = zeros(size(InitialBraL))';
    %TEST = 0;
    for BrasU = 1:length(AllBras)
        BrasU
        for KetsU = 1:length(AllBras)
            InitialBraU = AllBras(BrasU,:);
            InitialKetU = AllKets(:,KetsU);
            
    for deltaCountINT = 1:length(deltaINT)
        %TEST  =TEST +1;
        k1 = deltaINT(deltaCountINT,1);
        k2 = deltaINT(deltaCountINT,2);
        l1 = deltaINT(deltaCountINT,3);
        l2 = deltaINT(deltaCountINT,4);
        nk1 = deltaINT(deltaCountINT,9);
        nk2 = deltaINT(deltaCountINT,10);
        n11 = deltaINT(deltaCountINT,11);
        nl2 = deltaINT(deltaCountINT,12);
        mk1 = deltaINT(deltaCountINT,5);
        mk2 = deltaINT(deltaCountINT,6);
        ml1 = deltaINT(deltaCountINT,7);
        ml2 = deltaINT(deltaCountINT,8);
        OtherConst = deltaINT(deltaCountINT,16);
        
        l1Vect = l1Vect.*0;
        l2Vect = l2Vect.*0;
        k1Vect = k1Vect.*0;
        k2Vect = k2Vect.*0;
        
        l1Vect(l1,1) = -1;
        l2Vect(l2,1) = -1;
        k1Vect(k1,1) = 1;
        k2Vect(k2,1) = 1;
        
        l2Const = sqrt(InitialKetU(l2,1));
        l2trans  =InitialKetU + l2Vect;
        l1Const = sqrt(l2trans(l1,1));
        l1trans = l2trans + l1Vect;
        k2Const = sqrt(l1trans(k2,1)+1);
        k2trans = l1trans + k2Vect;
        k1Const = sqrt(k2trans(k1,1)+1);
        FinalKet = k2trans + k1Vect;
        
        if l2Const < 0;
            deltaCountINT
        end
          if l1Const < 0;
            deltaCountINT
        end
        
        
%                 TESTMAT(TEST,1) = l2Const;
%                 TESTMAT(TEST,2) = l1Const;
%                 TESTMAT(TEST,3) = k2Const;
%                 TESTMAT(TEST,4) = k1Const;
        
        if FinalKet == InitialBraU'
           % TEST = TEST + 1
            OpConst = l2Const*l1Const*k2Const*k1Const;
            %OtherConst = deltaINT(deltaCountINT,13);
           % TESTMAT(TEST,1) = OpConst;
            
            UMatNogConst(BrasU,KetsU) = UMatNogConst(BrasU,KetsU) + (OpConst.*OtherConst);
        end
    end
        end
    end
    
    %%%%% POT TERM
    if A > 0
    k1Vect = zeros(size(InitialBraL))';
    k2Vect = zeros(size(InitialBraL))';
    
    VMatNoA = zeros(length(AllKets));
    
    for BrasV = 1:length(AllBras)
        for KetsV = 1:length(AllBras)
            InitialBraV = AllBras(BrasV,:);
            InitialKetV = AllKets(:,KetsV);
    
            for deltaCountPOT = 1:length(deltaPOT)
                
                k1 = deltaPOT(deltaCountPOT,1);
                k2 = deltaPOT(deltaCountPOT,2);
                nk1 = deltaPOT(deltaCountPOT,3);
                nk2 = deltaPOT(deltaCountPOT,4);
                mk1 = deltaPOT(deltaCountPOT,5);
                mk2 = deltaPOT(deltaCountPOT,6);
                
                k1Vect = k1Vect.*0;
                k2Vect = k2Vect.*0;
                
                k1Vect(k1,1) = 1;
                k2Vect(k2,1) = -1;
                
                k2Const = sqrt(InitialKetV(k2,1));
                k2trans = InitialKetV + k2Vect;
                k1Const = sqrt(k2trans(k1,1)+1);
                FinalKetV = k2trans + k1Vect;
              
                
                
                
                if FinalKetV == InitialBraV'
                    
                    OpConst = k2Const*k1Const;
                    OtherConst = deltaPOT(deltaCountPOT,7);
                    VMatNoA(BrasV,KetsV) = VMatNoA(BrasV,KetsV) * (OtherConst*OpConst);
                end
            end
        end
    end
    else
    end

    
    n = 0;
for Omega = 0.6:0.0001:1.2;
    n=n+1;
    
    NAll = NMat;
    LAll = LMat*(1-Omega);
    UAll = UMatNogConst.*(g/(4*pi));
    UAlla=round(UAll,5);
    Total = NAll +UAll+(LAll);
    Eig = eig(Total);
    [V D] = eig(Total);
    
    
    OmegaIncreaseMatrix(1,n) = Omega;
    OmegaEnergyMatrix(:,n) = Eig;
    %%% Following line selects only the lowest 8 eigenstates
    OmegaMatRestricted = OmegaEnergyMatrix(1:9,:);
    
    
end

plot(OmegaIncreaseMatrix,OmegaEnergyMatrix);
xlabel('Omega')
ylabel('<E>')
ylim([8.2 8.9])
xlim([0.7 0.8])
toc
end
           
    



        
    
    
    
    
    
    
    
    
    
    
    

