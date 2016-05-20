clc; clear all;
tic
%%%%% DELETEMATRIX CHECK IF RECALCULATE IS ENABLED!!!!!!
DeleteMatrix = 1; %0 = no, 1 = yes
%%%%% DELETEMATRIX CHECK IF RECALCULATE IS ENABLED!!!!!!
N =6;
g = 1;
A = 0.00;
LLVal = 1; % LLL = 1, LL2 = 2
%ReCalculate = 0;% 0 = use matrix previously generated, 1 = do afresh

if LLVal == 1
    AllBras = csvread('StatesLLL.csv');
elseif LLVal ==2
    AllBras = csvread('StatesLLL_and_LL2.csv');
end
AllKets = AllBras';

KPos = [ 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19; %Matlab position
    0 0 1 0 1 0 1 0 1 0  1  0  1  0  1  0  1  0  1 ; %n value
    -1 0 0 1 1 2 2 3 3 4  4  5  5  6  6  7  7  8  8]; %Mt value

MomValue = [ -1 0 0 1 1 2 2 3 3 4 4 5 5 6 6 7 7 8 8]';
LengthAllBras = length(AllBras);
LengthVect = length(AllBras(1,:));

NMatFinal = (eye(LengthAllBras).*N);
csvwrite('NMatFinal.csv',NMatFinal)

for BraL = 1:length(AllBras)
    for KetL = 1:length(AllKets)
        InitialBraL = AllBras(BraL,:);
        InitialKetL = AllKets(:,KetL);
        if InitialBraL == InitialKetL'
            LMat(BraL,KetL) = dot(InitialKetL, MomValue);
        end
    end
end
csvwrite('LMatLLL.csv',LMat)

% Interaction Term, U
UMatNog = zeros(39);
k1Vect = zeros(19,1);
k2Vect = zeros(19,1);
l1Vect = zeros(19,1);
l2Vect = zeros(19,1);
%Udelta
if LLVal == 1;
    Ua = 0; Ub = 0; Uc = 0;
    for d1 = 2:2:18
        for d2 = 2:2:18
            for d3 = 2:2:18
                for d4 = 2:2:18
                    Ua = Ua + 1;
                    AllUdelta(Ua,1) = d1;
                    AllUdelta(Ua,2) = d2;
                    AllUdelta(Ua,3) = d3;
                    AllUdelta(Ua,4) = d4;
                    
                    if (AllUdelta(Ua,1) + AllUdelta(Ua,2)) == (AllUdelta(Ua,3) + AllUdelta(Ua,4))
                        Ub = Ub + 1;
                        UdeltaClean(Ub,1) = AllUdelta(Ua,1);
                        UdeltaClean(Ub,2) = AllUdelta(Ua,2);
                        UdeltaClean(Ub,3) = AllUdelta(Ua,3);
                        UdeltaClean(Ub,4) = AllUdelta(Ua,4);
                    end
                    
                end
            end
        end
    end
    %TESTNUMBER = 1;
    %for BrasU = TESTNUMBER;% 1:length(AllBras);
    %    for KetsU = TESTNUMBER;%1:length(AllBras);
            for BrasU = 1:length(AllBras);
        for KetsU = 1:length(AllBras);
            InitialBraU = AllBras(BrasU,:);
                InitialKetU = AllKets(:,KetsU);
                Uc = Uc.*0;
                for ScrollDelta = 1:length(UdeltaClean)
                                
                Uc = Uc + 1;
                
                k1 = UdeltaClean(Uc,1);
                k2 = UdeltaClean(Uc,2);
                l1 = UdeltaClean(Uc,3);
                l2 = UdeltaClean(Uc,4);
                
                nk1 = KPos(2,k1);
                nk2 = KPos(2,k2);
                nl1 = KPos(2,l1);
                nl2 = KPos(2,l2);
                
                mk1 = KPos(3,k1);
                mk2 = KPos(3,k2);
                ml1 = KPos(3,l1);
                ml2 = KPos(3,l2);
                
                k1Vect = k1Vect.*0;
                k2Vect = k2Vect.*0;
                l1Vect = l1Vect.*0;
                l2Vect = l2Vect.*0;
                
                k1Vect(k1,1) = 1;
                k2Vect(k2,1) = 1;
                l1Vect(l1,1) = -1;
                l2Vect(l2,1) = -1;
                
                l2Const = sqrt(InitialKetU(l2,1));
                l2trans = InitialKetU + l2Vect;
                l1Const = sqrt(l2trans(l2,1));
                l1trans = l2trans + l1Vect;
                k2Const = sqrt(l1trans(k2,1)+1);
                k2trans = l1trans + k2Vect;
                k1Const = sqrt(k2trans(k1,1)+1);
                FinalKet = k2trans +k2Vect;
                
                if FinalKet == InitialBraU'
                    PiTerm = factorial(nk1)*factorial(nk2)*factorial(nl1)*factorial(nl2)/...
                        factorial(nk1+abs(mk1))*factorial(nk2+abs(mk2))*factorial(nl1+abs(ml1))*factorial(nl2+abs(ml2));
                    
                    Mt = (abs(mk1)+abs(mk2)+abs(ml1)+abs(ml2));
                    MomTerm = 1/(2^(Mt/2));
                    
                    fun = @(x) exp(-x).*(x.^(Mt/2)).* (laguerreL(nk1,abs(mk1)).*(x/2)) .* (laguerreL(nk2,abs(mk2)).*(x/2)) .* (laguerreL(nl1,abs(ml1)).*(x/2)) .* (laguerreL(nl2,abs(ml2)).*(x/2));
                    
                    I2func = integral(fun,0,inf);
                    Const = l2Const*l1Const*k2Const*k1Const*PiTerm*MomTerm*I2func;
                    UMatNog(BrasU,KetsU) = UMatNog(BrasU,KetsU) + Const;
                end        
            end           
        end
    end
     
end

UMatFinal = (g/(4*pi))*UMatNog;
csvwrite('UMatFinalLLL.csv',UMatFinal);

if LLVal == 2;
    adasd
end

toc



if DeleteMatrix == 0
    delete('DeltaInteractionVariableLL.csv','LMatLLL.csv','NMatLLLFinal.csv','UMatFinalLLL.csv','VMatNoA_LLL.csv')
end