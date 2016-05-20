clc; clear all;
KPos = [1 2 3 4 5 6 7 8 9; %K
        0 1 2 3 4 5 6 7 8];%Mt

    
    aa=0;
for k1 = 1:9
    for k2 = 1:9
        for k3 = 1:9
            for k4 = 1:9
                aa = aa+1;
                All(aa,1) = k1;
                All(aa,2) = k2;
                All(aa,3) = k3;
                All(aa,4) = k4;
            end
        end
    end
end
bb=0;
for CountIntAll = 1:length(All)
    k1 = All(CountIntAll,1);
    k2 = All(CountIntAll,2);
    k3 = All(CountIntAll,3);
    k4 = All(CountIntAll,4);
    if (KPos(2,k1) + KPos(2,k2)) == (KPos(2,k3) + KPos(2,k4))
        bb=bb+1;
        IntDelta(bb,:) = All(CountIntAll,:);    
    end
end
csvwrite('IntDelta.csv',IntDelta);

cc=0;
for Pk1 = 1:9;
    for Pk2 = 1:9;
        cc=cc+1;
        PAll(cc,1)=Pk1;
        PAll(cc,2)=Pk2;
    end
end
dd=0;
for CountPotAllP = 1:length(PAll)
    Pk1 = PAll(CountPotAllP,1);
    Pk2 = PAll(CountPotAllP,2);
    if Pk1 == (Pk2 + 2)
        dd=dd+1;
        PotDeltaP(dd,1) = Pk1;
        PotDeltaP(dd,2) = Pk2;
    end
end
ee=0;
for CountPotAllM = 1:length(PAll)
    Pk1 = PAll(CountPotAllM,1);
    Pk2 = PAll(CountPotAllM,2);
    if Pk1 == (Pk2 - 2)
        ee=ee+1;
        PotDeltaM(ee,1) = Pk1;
        PotDeltaM(ee,2) = Pk2;
    end
end
PotDelta = [PotDeltaP; PotDeltaM];
csvwrite('PotDelta.csv',PotDelta);

        
    
