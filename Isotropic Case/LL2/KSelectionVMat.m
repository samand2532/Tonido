function [ AllowedKVal ] = KSelectionVMat( InitialKet )
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here

KPos = [ 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22; %Matlab position
         0 0 0 0 0 0 0 0 0 0  0  1  1  1  1  1  1  1  1  1  1  1 ; %n value
        -1 0 1 2 3 4 5 6 7 8  9 -1  0  1  2  3  4  5  6  7  8  9]; %Mt va
    
KPos = [ 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22; %Matlab position
         0 0 0 0 0 0 0 0 0 0  0  1  1  1  1  1  1  1  1  1  1  1 ; %n value
        -1 0 1 2 3 4 5 6 7 8  9 -1  0  1  2  3  4  5  6  7  8  9]; %Mt va
    
A = InitialKet;
aa=0;bb=0;cc=0;dd=0;ee=0;ff=0;gg=0;hh=0;kk=0;ll=0;
subVect = zeros(1,length(A));
%%%%%annihilation part

for CountA = 1:length(A)
    subVect = subVect.*0;
    subVect(1,CountA) = -1;
    sub = A + subVect;
    aa = aa + 1;
    MatA(aa,:) = sub;
    deltaP1(aa,1) = CountA;
end
[x,y] = size(MatA);    
for CountB = 1:x
    if any(MatA(CountB,:) < 0)
    else bb = bb +1;
        MatB(bb,:) = MatA(CountB,:);
        deltaP2(bb,:) = deltaP1(CountB,:);
    end
end
%%%%%%%%%%%CreationPart

addVect = zeros(1,length(A));
[x,y] = size(MatB);
for CountC =1:x
    initialVect = MatB(CountC,:);
    k2 = deltaP2(CountC,1);
    for CountD = 1:length(initialVect)
        addVect = addVect.*0;
        addVect(1,CountD) = 1;
        add = initialVect + addVect;
        cc = cc + 1;
        MatC(cc,:) = add;
        delta(cc,2) = k2;
        delta(cc,1) = CountD;
    end
end
[x,y] = size(delta);
for CountD = 1:x
    k1 = delta(CountD,1);
    k2 = delta(CountD,2);
    
    delta(CountD,3) = KPos(2,k1);
    delta(CountD,4) = KPos(2,k2);
    
    delta(CountD,5) = KPos(3,k1);
    delta(CountD,6) = KPos(3,k2);
end
[x,y] = size(delta);
for CountE = 1:x
    if delta(CountE,5) == (delta(CountE,6) + 2) || delta(CountE,5) == (delta(CountE,6) - 2)
        dd = dd + 1;
        deltaMt(dd,:) = delta(CountE,:);
    end
end

[x,y] = size(deltaMt);
for CountF = 1:x
    if (delta(CountF,3) + delta(CountF,4)) < 2
        ee = ee + 1;
        deltaNMt(ee,:) = deltaMt(CountF,:);
    end
end

[x,y] = size(deltaNMt);

for CountG = 1:x
    k2 = deltaNMt(CountG,2);
    k1 = deltaNMt(CountG,1);
    subVect = subVect.*0;
    subVect(1,k2) = -1;
    addVect = addVect.*0;
    addVect(1,k1) = 1;
    Final = A + subVect + addVect;
    
    if (dot(Final,KPos(2,:)) < 2) && (Final(1,1) < 2) && (Final(1,12) == 0) && sum(Final(1,1)+any(Final(1,12:22)))<2 ;
        ff = ff +1;
        AllowedKVal(ff,:) = deltaNMt(CountG,:);
    end
    
     
end





end

