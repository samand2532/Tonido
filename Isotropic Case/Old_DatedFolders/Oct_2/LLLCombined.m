tic
clear all; clc;
DeleteMatrix = 0; % 1 = yes
N = 6; g =1; LLVal = 1;% For the deltainteractionfunction, LLL = 1, LL2incLLL = 2
A = 0.03; 
MomValue = [ -1 0 0 1 1 2 2 3 3 4 4 5 5 6 6 7 7 8 8]'; 
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

%%%% INteraction Term
delta = csvread('LagIntLLL.csv');
k1Vect = zeros(19,1);
k2Vect = zeros(19,1);
l1Vect = zeros(19,1);
l2Vect = zeros(19,1);

IntMat = zeros(39);

for BrasU = 1:length(AllBras)
    for KetsU = 1:length(AllKets)
        InitialBraU = AllBras(BrasU,:);
        InitialKetU = AllKets(:,KetsU);
        
        for CountDelta = 1:length(delta)
            
            k1Vect = k1Vect.*0;
            k2Vect = k2Vect.*0;
            l1Vect = l1Vect.*0;
            l2Vect = l2Vect.*0;
            
            nk1 = delta(CountDelta,1);
            nk2 = delta(CountDelta,2);
            nl1 = delta(CountDelta,3);
            nl2 = delta(CountDelta,4);
            mk1 = delta(CountDelta,5);
            mk2 = delta(CountDelta,6);
            ml1 = delta(CountDelta,7);
            ml2 = delta(CountDelta,8);
            
            k1Vect(mk1,1) = 1;
            k2Vect(mk2,1) = 1;
            l1Vect(ml1,1) = -1;
            l2Vect(ml2,1) = -1;
            
            l2Const = sqrt(InitialKetU(ml2,1));
            l2trans = InitialKetU + l2Vect;
            l1Const = sqrt(l2trans(ml1,1));
            l1trans = l2trans + l1Vect;
            k2Const = sqrt(l1trans(mk2,1)+1);
            k2trans = l1trans + k2Vect;
            k1Const = sqrt(l2trans(ml1,1)+1);
            FinalKet = k2trans + k1Vect;
            
            if FinalKet == InitialBraU'
                OpConst = l2Const*l1Const*k2Const*k1Const
                OtherVariables = delta(CountDelta,9)
            end
            IntMat(BraU,KetsU) = IntMat(BraU,KetsU) + (OpConst*OtherVariables);
        
        end
 
    end
end

