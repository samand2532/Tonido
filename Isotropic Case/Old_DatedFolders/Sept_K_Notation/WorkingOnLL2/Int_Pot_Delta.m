
clc;clear all;

LengthVect =19;
MomValue = [ -1 0 0 1 1 2 2 3 3 4 4 5 5 6 6 7 7 8 8];  

KPos = [ 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19; %Matlab position
         0 0 1 0 1 0 1 0 1 0  1  0  1  0  1  0  1  0  1 ; %n value
        -1 0 0 1 1 2 2 3 3 4  4  5  5  6  6  7  7  8  8]; %Mt value

aa=0;    
StorageMatrixAllK = zeros(LengthVect^4,4);
for k1 = 1:LengthVect;
    for k2 = 1:LengthVect;
        for l1 = 1:LengthVect;
            for l2 = 1:LengthVect;
                aa = aa+1;
                StorageMatrixAllK(aa,1) = k1;
                StorageMatrixAllK(aa,2) = k2;
                StorageMatrixAllK(aa,3) = l1;
                StorageMatrixAllK(aa,4) = l2;
            end
        end
    end
end
bb=0;
%%%% InteractionTermDelta
%%%% Section balances momentum totals
for Countm = 1:length(StorageMatrixAllK);
    k1 = StorageMatrixAllK(Countm,1);
    k2 = StorageMatrixAllK(Countm,2);
    l1 = StorageMatrixAllK(Countm,3);
    l2 = StorageMatrixAllK(Countm,4);
if (KPos(3,k1) + KPos(3,k2)) == (KPos(3,l1) + KPos(3,l2));
    bb=bb+1;
    IntDeltaMomBalanced(bb,:) = StorageMatrixAllK(Countm,:);
end
end
%%%% section removes multiple n's
cc=0;
for Countn = 1:length(IntDeltaMomBalanced);
    k1 = IntDeltaMomBalanced(Countn,1);
    k2 = IntDeltaMomBalanced(Countn,2);
    l1 = IntDeltaMomBalanced(Countn,3);
    l2 = IntDeltaMomBalanced(Countn,4);
    if (KPos(2,k1) + KPos(2,k2)) + (KPos(2,l1) + KPos(2,l2)) > 1
    else cc = cc+1;
        IntDeltaMomNBal(cc,:) = IntDeltaMomBalanced(Countn,:);
    end
end
%%%%% section removes multiple Mt neg's
dd= 0;
removeMat = (IntDeltaMomNBal == 1);
removeMatsum = sum(removeMat,2);
for CountRMS = 1:length(removeMatsum);
    if removeMatsum(CountRMS,1) < 2;
        dd=dd+1;
        IntDelta(dd,:) = IntDeltaMomNBal(CountRMS,:);
    end
end
%%%% IntDelta is all allowwed K values, for LL2 basis (includes LLL)
csvwrite('IntDelta.csv',IntDelta);


%%%%%%%% Potential Term Delta
ee=0;

for Pk1 = 1:19;
    for Pk2 = 1:19;
        ee=ee+1;
        PInitialMat(ee,1) = Pk1;
        PInitialMat(ee,2) = Pk2;
    end
end
%%%% Section removes multiple n's
ff=0;
for PCountn = 1:length(PInitialMat);
    Pk1 = PInitialMat(PCountn,1);
    Pk2 = PInitialMat(PCountn,2);
    if (KPos(2,Pk1) + KPos(2,Pk2)) > 1;
    else ff = ff+1;
        PDeltaNBal(ff,:) = PInitialMat(PCountn,:);
    end
end
PremoveMat = (PDeltaNBal == 1);
PremoveMatSum = sum(PremoveMat,2);
gg = 0;
for CountPRMS = 1:length(PremoveMatSum);
    if PremoveMatSum(CountPRMS,1) < 2;
        gg=gg+1;
        PotDelta(gg,:) = PDeltaNBal(CountPRMS,:);
    end
end

hh=0;
for CountPotDeltaP = 1:length(PotDelta)
    Pk1 = PotDelta(CountPotDeltaP,1);
    Pk2 = PotDelta(CountPotDeltaP,2);
    if KPos(3,Pk1) == (KPos(3,Pk2) + 2)
        hh = hh + 1;
        PotDeltaP(hh,1) = PotDelta(CountPotDeltaP,1);
        PotDeltaP(hh,2) = PotDelta(CountPotDeltaP,2);
    end
end
ii = 0;
for CountPotDeltaM = 1:length(PotDelta)
    Pk1 = PotDelta(CountPotDeltaM,1);
    Pk2 = PotDelta(CountPotDeltaM,2);
    if KPos(3,Pk1) == (KPos(3,Pk2) - 2)
        ii = ii + 1;
        PotDeltaM(ii,1) = PotDelta(CountPotDeltaM,1);
        PotDeltaM(ii,2) = PotDelta(CountPotDeltaM,2);
    end
end

PotDeltaPM = [PotDeltaP; PotDeltaM];
csvwrite('PotDeltaPM.csv',PotDeltaPM);
        

        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        