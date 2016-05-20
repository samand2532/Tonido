function [] = AnisoFunc()
%ANISOFUNC Summary of this function goes here
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
    
%%%%%% SEction makes +/-2 

CountC = 0;
CountE = 0;
for CountD = 1:length(AllKMat)
    k1 = AllKMat(CountD,1);
        k2 = AllKMat(CountD,2);
        l1 = AllKMat(CountD,3);
        l2 = AllKMat(CountD,4);
        
        mk1=KPos(2,k1);mk2=KPos(2,k2);ml1=KPos(2,l1);ml2=KPos(2,l2);
        
        if (mk1+mk2) == ((ml1+ml2)+2)
            CountC = CountC + 1;
            MatAa(CountC,:) = AllKMat(CountD,:);
        end
        if (mk1+mk2) == ((ml1+ml2)-2)
            CountE = CountE + 1;
            MatAb(CountE,:) = AllKMat(CountD,:);
        end
end
MatA = [MatAa;MatAb];
        %%%%% Stop mk or ml > 8
 CountF = 0;
 for CountG = 1:length(MatA)
        k1 = MatA(CountG,1);
        k2 = MatA(CountG,2);
        l1 = MatA(CountG,3);
        l2 = MatA(CountG,4);
        
        mk1=KPos(2,k1);mk2=KPos(2,k2);ml1=KPos(2,l1);ml2=KPos(2,l2);
        
        if ((mk1 + mk2) <9) && ((ml1+ml2) < 9)
            CountF = CountF + 1;
            MatB(CountF,:) = MatA(CountG,:);
        end
 end
 
 nMat = MatB(:,5:8) == 1;
sumnMat = sum(nMat,2);

CountH = 0;
for CountI = 1:length(MatB)
    if sumnMat(CountI) < 2
        CountH = CountH + 1;
        MatC(CountH,:) = MatB(CountI,:);
    end
end

mMat = MatC(:,9:12) == -1;
summMat = sum(mMat,2);

CountJ = 0;
for CountK = 1:length(MatC)
    if summMat(CountK) < 2
        CountJ = CountJ + 1;
        MatD(CountJ,:) = MatC(CountK,:);
    end
end

nMatA = MatD(:,5:8) == 1;
mMatA = MatD(:,9:12) == -1;
nmMatA = [nMatA mMatA];
sumnmMatA = sum(nmMatA,2);

CountL = 0;
for CountM = 1:length(MatD)
    if sumnmMatA(CountM) < 2
        CountL = CountL + 1;
        MatE(CountL,:) = MatD(CountM,:);
    end
end
    %%%% 1  2  3  4  5   6   7   8    9  10   11  12    13     14
    %%%% k1|k2|l1|l2|nk1|nk2|nl1|nl2|mk1|mk2|ml1|ml2|SqrtTerm|I1
for CountN = 1:length(MatE)
    CountN
    k1 = MatE(CountN,1);
    k2 = MatE(CountN,2);
    
    mk1 = KPos(2,k1); mk2 = KPos(2,k2);
    nk1 = KPos(3,k1); nk2 = KPos(3,k2);
    
    SqrtTerm = sqrt(1/(factorial(nk1+abs(mk1))*(factorial(nk2+abs(mk2)))));
   MatE(CountN,13) = SqrtTerm;
   Mt = (abs(mk1) + abs(mk2) + 2)/2;
   fun = @(x) exp(-x).*(x.^(Mt)).*(laguerreL(nk1,abs(mk1),(x/2))).*(laguerreL(nk2,abs(mk2),(x/2)));
   I1func = integral(fun,0,inf);
   
   MatE(CountN,14) = I1func;
   MatE(CountN,15) = SqrtTerm*I1func;
end

AMatNoConst = zeros(length(AllBras));

            k2Vect = zeros(19,1);
            k1Vect = zeros(19,1);

for BrasA = 1:length(AllBras)
    InitialBra = AllBras(BrasA,:)
    for KetsA = 1:length(AllBras)
        InitialKet = AllKets(:,KetsA);
        
        for Delta = 1:length(MatE)
            k2Vect = k2Vect.*0;
            k1Vect = k1Vect.*0;
            
            k1 = MatE(Delta,1);
            k2 = MatE(Delta,2);
            
            k1Vect(k1,1) = 1;
            k2Vect(k2,1) = -1;

            k2Const = sqrt(InitialKet(k2,1));
            k2trans = InitialKet + k2Vect;
            k1Const = sqrt((k2trans(k1,1))+1);
            FinalKet = k2trans + k1Vect;
            
            if InitialBra == FinalKet'
                OpConst = k2Const*k1Const;
                OtherConst = MatE(Delta,15);
                Tot = OpConst*OtherConst;
                AMatNoConst(BrasA,KetsA) = AMatNoConst(BrasA,KetsA) + Tot;
                
            end
            
            
            
        end
    end
end




csvwrite('AnisotermNoConst.csv',AMatNoConst);

end

