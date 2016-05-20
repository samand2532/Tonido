clc; clear all; close all;
tic
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
%%%% Following limit deltaMt < 9
CountB = 0;
for CountA = 1:length(FullMat)
    
    if (FullMat(CountA,9)+ FullMat(CountA,10)) < 9
        CountB = CountB + 1;
        FullMatA(CountB,:) = FullMat(CountA,:);
    end
end

%%%%% section makes mk = ml + 2

summk = sum(FullMatA(:,(9:10)),2);
summl = sum(FullMatA(:,(11:12)),2);

CountAa = 0;
for CountAb = 1:length(FullMatA)
    if summk(CountAb) == (summl(CountAb) + 2)
        CountAa = CountAa+1;
        FullMatAa(CountAa,:) = FullMatA(CountAb,:);
    elseif summk(CountAb) == (summl(CountAb) - 2)
        CountAa = CountAa + 1;
        FullMatAa(CountAa,:) = FullMatA(CountAb,:);
    end
end

CountAc = 0;
for CountAd = 1:length(FullMatAa)
    if (FullMatAa(CountAd,9) + FullMatAa(CountAd,10)) < 9 &&...
            (FullMatAa(CountAd,11) + FullMatAa(CountAd,12)) < 9
        
        CountAc = CountAc + 1;
        FullMatB(CountAc,:) = FullMatAa(CountAd,:);
    end
end


%%%%% delets multile n's
deleteN = FullMatB(:,5:8) == 1;
sumDelN = sum(deleteN,2);

CountE = 0;
for CountF = 1:length(FullMatB)
    if sumDelN(CountF) < 2
        CountE = CountE + 1;
        FullMatC(CountE,:) = FullMatB(CountF,:);
    end
end
      
%%%%% delets multile mt=-1s
deleteM = FullMatC(:,9:12) == -1;
sumDelM = sum(deleteM,2);

CountG = 0;
for CountH = 1:length(FullMatC)
    if sumDelM(CountH) < 2
    
        CountG = CountG + 1;
        FullMatD(CountG,:) = FullMatC(CountH,:);
    end
end        

%%%%% stops mt=-1 and n=1 occuring in same line

sumN = sum(FullMatD(:,5:8),2);
delMt = FullMatD(:,9:12) == -1;
sumM = sum(delMt,2);

CountI = 0;
for CountJ = 1:length(FullMatD)
    if (sumM(CountJ) + sumN(CountJ)) < 2
        CountI = CountI + 1;
        AnisoDelta(CountI,:) = FullMatD(CountJ,:);
    end
end

        
toc
        
        
        
        
        
