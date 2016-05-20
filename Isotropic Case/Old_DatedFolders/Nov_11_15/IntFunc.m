function [] = IntFunc()
%INTFUNC Summary of this function goes here
%   Detailed explanation goes here
AllBras = csvread('NewLL2.csv');
AllKets = AllBras';

KPos = [ 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19;%%K number
    -1 0 0 1 1 2 2 3 3 4  4  5  5  6  6  7  7  8  8;%%%Mt
    0 0 1 0 1 0 1 0 1 0  1  0  1  0  1  0  1  0  1];%%% n value

AllKMat = zeros(130321,4);
  
    CountA = 0;
    for k1 = 1:19
        for k2 = 1:19
            for k3 = 1:19
                for k4 = 1:19
                    CountA = CountA + 1;
                    AllKMat(CountA,1) = k1;
                    AllKMat(CountA,2) = k2;
                    AllKMat(CountA,3) = k3;
                    AllKMat(CountA,4) = k4;
                end
            end
        end
    end
    %%%% 1  2  3  4  5   6   7   8    9  10   11  12
    %%%% k1|k2|l1|l2|nk1|nk2|nl1|nl2|mk1|mk2|ml1|ml2
    
    for CountB = 1:length(AllKMat)
        k1 = AllKMat(CountB,1);
        k2 = AllKMat(CountB,2);
        l1 = AllKMat(CountB,3);
        l2 = AllKMat(CountB,4);
        
        AllKMat(CountB,5) = KPos(3,k1);
        AllKMat(CountB,6) = KPos(3,k2);
        AllKMat(CountB,7) = KPos(3,l1);
        AllKMat(CountB,8) = KPos(3,l2);
        AllKMat(CountB,9) = KPos(2,k1);
        AllKMat(CountB,10) = KPos(2,k2);
        AllKMat(CountB,11) = KPos(2,l1);
        AllKMat(CountB,12) = KPos(2,l2);
    end
    
    %%%%%% Section balances Mk and Ml
    CountBb =0;
    for CountBa = 1:length(AllKMat)
        k1 = AllKMat(CountBa,1);
        k2 = AllKMat(CountBa,2);
        l1 = AllKMat(CountBa,3);
        l2 = AllKMat(CountBa,4);
        
        mk1=KPos(2,k1);mk2=KPos(2,k2);ml1=KPos(2,l1);ml2=KPos(2,l2);
        
        if ((mk1+mk2) == (ml1+ml2))
            CountBb = CountBb + 1;
            AllKMata(CountBb,:) = AllKMat(CountBa,:);
        end
        
    end
    %%%% Section remove Mt=-1's and n=1's and Mt and n occuring together
    
    
    mMat = AllKMata(:,9:12) == -1;
    summMat = sum(mMat,2);
    
    CountC = 0;
    for CountD = 1:length(AllKMata)
        if (summMat(CountD) < 2)
            CountC = CountC + 1;
            MatA(CountC,:) = AllKMata(CountD,:);
        end
    end
    
    nMat = MatA(:,5:8) == 1;
    sumnMat = sum(nMat,2);
    
    CountE = 0;
    for CountF = 1:length(MatA)
        if (sumnMat(CountF) > 1)
        else
            CountE = CountE + 1;
            MatB(CountE,:) = MatA(CountF,:);
        end
    end
    
    mMatA = MatB(:,9:12) == -1;
    nMatA = MatB(:,5:8) == 1;
    CombMatA = [mMatA nMatA];
    sumCombMat = sum(CombMatA,2);
    
    CountG = 0;
    for CountH = 1:length(MatB)
        if sumCombMat(CountH) > 1
        else CountG = CountG + 1;
            MatC(CountG,:) = MatB(CountH,:);
        end
    end
    
    CountI = 0;
    for CountJ = 1:length(MatC)
        k1 = MatC(CountJ,1); k2 = MatC(CountJ,2);
        mk1=KPos(2,k1);mk2=KPos(2,k2);
        if (mk1 + mk2) < 9
            CountI = CountI + 1;
            Delta(CountI,:) = MatC(CountJ,:);
        end
    end
    
    %%%% 1  2  3  4  5   6   7   8    9  10   11  12|   13     14    15
    %%%% k1|k2|l1|l2|nk1|nk2|nl1|nl2|mk1|mk2|ml1|ml2|MomTerm|PiTerm|I2
    
    for CountK = 1:length(Delta)
        
        k1 = Delta(CountK,1);
        k2 = Delta(CountK,2);
        l1 = Delta(CountK,3);
        l2 = Delta(CountK,4);
        
        mk1=KPos(2,k1);mk2=KPos(2,k2);ml1=KPos(2,l1);ml2=KPos(2,l2);
        nk1=KPos(3,k1);nk2=KPos(3,k2);nl1=KPos(3,l1);nl2=KPos(3,l2);
        Mt = abs(mk1)+abs(mk2)+abs(ml1)+abs(ml2);
        
        fun = @(x) exp(-x).*(x.^(Mt)).*(laguerreL(nk1,abs(mk1),(x/2))).*(laguerreL(nk2,abs(mk2),(x/2))).*(laguerreL(nl1,abs(ml1),(x/2))).*(laguerreL(nl2,abs(ml2),(x/2)));
        I2func = integral(fun,0,inf);
        
        Delta(CountK,13) = 1/(2^(((abs(Delta(CountK,9)))+(abs(Delta(CountK,10)))+(abs(Delta(CountK,11)))+(abs(Delta(CountK,12))))/2));
        Delta(CountK,14) = sqrt(1/((factorial(nk1+abs(mk1)))*(factorial(nk2+abs(mk2)))*(factorial(nl1+abs(ml1)))*(factorial(nl2+abs(ml2)))));
        Delta(CountK,15) = I2func;
        Delta(CountK,16) = (Delta(CountK,13))*(Delta(CountK,14))*(Delta(CountK,15));
        
    end



%%%% Crossing Terms
UMatNoConst = zeros(length(AllBras));

            l2Vect = zeros(19,1);
            l1Vect = zeros(19,1);
            k2Vect = zeros(19,1);
            k1Vect = zeros(19,1);

for BrasU = 1:length(AllBras)
    BrasU
    for KetsU = 1:length(AllKets)
        %KetsU
        InitialBra = AllBras(BrasU,:);
        InitialKet = AllKets(:,KetsU);
        
        for CountIntDelta = 1:length(Delta)
            l2Vect = l2Vect.*0;
            l1Vect = l1Vect.*0;
            k2Vect = k2Vect.*0;
            k1Vect = k1Vect.*0;
            
            k1 = Delta(CountIntDelta,1);
            k2 = Delta(CountIntDelta,2);
            l1 = Delta(CountIntDelta,3);
            l2 = Delta(CountIntDelta,4);
            
            l2Vect(l2,1) = -1;
            l1Vect(l1,1) = -1;
            k1Vect(k1,1) = 1;
            k2Vect(k2,1) = 1;
            
            l2Const = sqrt(InitialKet(l2,1));
            l2trans = InitialKet + l2Vect;
            l1Const = sqrt(l2trans(l1,1));
            l1trans = l2trans + l1Vect;
            k2Const = sqrt(l1trans(k2,1)+1);
            k2trans = l1trans + k2Vect;
            k1Const = sqrt(k2trans(k1,1)+1);
            FinalKet = k2trans + k1Vect;
            
            if FinalKet == InitialBra'
                
                OpConst = l2Const*l1Const*k2Const*k1Const;
                OtherConst = Delta(CountIntDelta,16);
                TotalConst = OpConst*OtherConst;
                
                UMatNoConst(BrasU,KetsU) = UMatNoConst(BrasU,KetsU) + TotalConst;
            end
        end
    end

end
csvwrite('UMatNoConst.csv',UMatNoConst);
end

   
    
    



