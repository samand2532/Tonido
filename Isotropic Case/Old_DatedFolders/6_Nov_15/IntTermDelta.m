clc; clear all; close all;
tic

%%%%% SECTION AT END FOR CONSTS
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
  %%% 1  2  3  4  5   6    7  8   9  10  11  12
%%%% k1|k2|l1|l3|nk1|nk2|nl1|nl2|mk1|mk2|ml1|ml2

FullMat(:,:) = AllKMat(:,:);
FullMat(:,5) = KPos(3,FullMat(:,1));
FullMat(:,6) = KPos(3,FullMat(:,2));
FullMat(:,7) = KPos(3,FullMat(:,3));
FullMat(:,8) = KPos(3,FullMat(:,4));
FullMat(:,9) = KPos(2,FullMat(:,1));
FullMat(:,10) = KPos(2,FullMat(:,2));
FullMat(:,11) = KPos(2,FullMat(:,3));
FullMat(:,12) = KPos(2,FullMat(:,4));

%%%% remove multiple n's

CountB = 0;
for CountC = 1:length(FullMat)
    if sum(FullMat(CountC,5:8),2) > 1
    else CountB = CountB +1;
        Matrestrictn(CountB,:) = FullMat(CountC,:);
    end
end

%%%% remove multiple Mt=-1's
CountE = 0;
TempMat = Matrestrictn(:,9:12) == -1;
for CountD = 1:length(Matrestrictn)
    if sum(TempMat(CountD,:),2) > 1
    else CountE = CountE+1;
        MatRestrictMN(CountE,:) = Matrestrictn(CountD,:);
    end
end

%%%% remove Mt=-1 and n appearing on same line

sumTempMatn = sum((MatRestrictMN(:,5:8) == 1),2);
sumTempMatm = sum((MatRestrictMN(:,9:12) == -1),2);
VarA = sum([sumTempMatn sumTempMatm],2);

CountF = 0;
for CountG = 1:length(MatRestrictMN)
    if VarA(CountG) > 1
    else CountF = CountF + 1;
        MatRestrict(CountF,:) = MatRestrictMN(CountG,:);
    end
end

CountH = 0;
for CountI = 1:length(MatRestrict)
    if (KPos(2,MatRestrict(CountI,1))) + (KPos(2,MatRestrict(CountI,2))) < 9
        CountH = CountH +1;
        MatRestricted1(CountH,:) = MatRestrict(CountI,:);
    end
end

CountJ =0;
for CountK = 1:length(MatRestricted1)
    if ((KPos(2,MatRestricted1(CountK,1)))+(KPos(2,MatRestricted1(CountK,2)))) ...
            == ((KPos(2,MatRestricted1(CountK,3)))+(KPos(2,MatRestricted1(CountK,4))))
        CountJ = CountJ + 1;
        IntDelta(CountJ,:) = MatRestricted1(CountK,:);
    end
end

%%%% THIS SECTION IS FOR THE MOM,PiandI2CONSTS
%%%%% 1  2  3  4  5   6    7  8   9  10  11  12
%%%% k1|k2|l1|l3|nk1|nk2|nl1|nl2|mk1|mk2|ml1|ml2
for CountL = 1:length(IntDelta)
    nk1 = IntDelta(CountL,5);nk2 = IntDelta(CountL,6);nl1 = IntDelta(CountL,7);nl2 = IntDelta(CountL,8);
    mk1 = IntDelta(CountL,9);mk2 = IntDelta(CountL,10);ml1 = IntDelta(CountL,11);ml2 = IntDelta(CountL,12);
    Mt = abs(IntDelta(CountL,9))+abs(IntDelta(CountL,10))+abs(IntDelta(CountL,11))+abs(IntDelta(CountL,12));
    MomTerm = 1/(2^(Mt/2));
    PiTerm = 1/(factorial(nk1 +abs(mk1))*factorial(nk2 +abs(mk2))*factorial(nl1 +abs(ml1))*factorial(nl2 +abs(ml2)));
    fun = @(x) exp(-x).*(x.^(Mt/2)).* (laguerreL(nk1,abs(mk1)).*(x/2)) .* (laguerreL(nk2,abs(mk2)).*(x/2)) .* (laguerreL(nl1,abs(ml1)).*(x/2)) .* (laguerreL(nl2,abs(ml2)).*(x/2));
    I2func = integral(fun,0,inf);
    IntDelta(CountL,13) = MomTerm*PiTerm*I2func;
    
end
csvwrite('IntDeltainConst.csv', IntDelta);
   toc     
    
    

    

