function [ deltaNearFinal ] = KSelection( InitialKet )
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here
KPos = [ 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22; %Matlab position
         0 0 0 0 0 0 0 0 0 0  0  1  1  1  1  1  1  1  1  1  1  1 ; %n value
        -1 0 1 2 3 4 5 6 7 8  9 -1  0  1  2  3  4  5  6  7  8  9]; %Mt va
    
A = InitialKet;
aa=0;bb=0;cc=0;dd=0;ee=0;ff=0;gg=0;hh=0;kk=0;ll=0;mm=0;nn=0;pp=0;qq=0;rr=0;
firstsubVect = zeros(1,length(A));
twosubVect = zeros(1,length(A));
%%%%%%%%%% Subtraction Section
for CountA = 1:length(A)
    firstsubVect = firstsubVect.*0;
    firstsubVect(1,CountA) = -1;
    firstSub = A+firstsubVect;
    for CountB = 1:length(A)
        twosubVect = twosubVect.*0;
        twosubVect(1,CountB) = -1;
        twosub = firstSub+twosubVect;
        aa = aa + 1;
        MatA(aa,:) = twosub;
        deltaP1(aa,1) = CountA;
        deltaP1(aa,2) = CountB;
    end
end

for CountC = 1:length(MatA)
    if any(MatA(CountC,:) < 0)
    else bb = bb + 1;
        MatB(bb,:) = MatA(CountC,:);
        deltaP2(bb,:) = deltaP1(CountC,:);
    end
end

%%%%%%%%%%%%%% AdditionSection
firstaddVect = zeros(1,length(A));
twoaddVect = zeros(1,length(A));        
[x y] = size(MatB);
for CountD = 1:x
    initialVect = MatB(CountD,:);
    k1 = deltaP2(CountD,1);
    k2 = deltaP2(CountD,2);
    for CountE = 1:length(initialVect)
        firstaddVect = firstaddVect.*0;
        firstaddVect(1,CountE) = 1;
        firstadd = initialVect+firstaddVect;
        for CountF = 1:length(initialVect)
            twoaddVect = twoaddVect.*0;
            twoaddVect(1,CountF) = 1;
            twoadd = firstadd + twoaddVect;
            cc = cc + 1;
            MatC(cc,:) = twoadd;
            delta(cc,1) = k1;
            delta(cc,2) = k2;
            delta(cc,3) = CountE;
            delta(cc,4) = CountF;
        end
    end
end

for CountG = 1:length(delta)
    k1 = delta(CountG,1);
    k2 = delta(CountG,2);
    l1 = delta(CountG,3);
    l2 = delta(CountG,4);
    
    delta(CountG,5)  = KPos(2,k1);
    delta(CountG,6) = KPos(2,k2);
    delta(CountG,7) = KPos(2,l1);
    delta(CountG,8) = KPos(2,l2);
    
    delta(CountG,9)  = KPos(3,k1);
    delta(CountG,10) = KPos(3,k2);
    delta(CountG,11) = KPos(3,l1);
    delta(CountG,12) = KPos(3,l2); 
end
    
for CountH = 1:length(delta)
    if (delta(CountH,9) + delta(CountH,10)) == (delta(CountH,11) + delta(CountH,12))
        hh = hh + 1;
        deltaMt(hh,:) = delta(CountH,:);
    end
end

[x y] = size(deltaMt);
for CountK = 1:x
    if ((deltaMt(CountK,5)+deltaMt(CountK,6)) < 2) && ((deltaMt(CountK,7)+deltaMt(CountK,8)) < 2)
        kk = kk + 1;
        deltaNMt(kk,:) = deltaMt(CountK,:);        
    end
end

[x y] = size(deltaNMt);
for CountL = 1:x
    if ((deltaNMt(CountL,9) == -1) && (deltaNMt(CountL,10) == -1)) ||...
            ((deltaNMt(CountL,11) == -1) && (deltaNMt(CountL,12) == -1))
    else ll = ll + 1;
        deltaNMtOne(ll,:) = deltaNMt(CountL,:);
    end
end

[x y] = size(deltaNMtOne);
for CountM = 1:x
    if (((deltaNMtOne(CountM,5) == 1) || (deltaNMtOne(CountM,6) == 1)) && ((deltaNMtOne(CountM,9) == -1) || (deltaNMtOne(CountM,10) == -1)))...
            || (((deltaNMtOne(CountM,7) == 1) || (deltaNMtOne(CountM,8) == 1)) && ((deltaNMtOne(CountM,11) == -1) || (deltaNMtOne(CountM,12) == -1)))
    else mm = mm +1;
        deltaNMtOneNoMtLim(mm,:) = deltaNMtOne(CountM,:);
    end    
end

[x y] = size(deltaNMtOneNoMtLim);
for CountN = 1:x
    if (deltaNMtOneNoMtLim(CountN,11) + deltaNMtOneNoMtLim(CountN,12)) < 10
        nn = nn + 1;
        deltaNearFinal(nn,:) = deltaNMtOneNoMtLim(CountN,:);
    end
end
%%%%% This section stops n values being created(or mt = -1) when the
%%%%% original state already has them. 
[x y] = size(deltaNearFinal);
BalanceVect = zeros(1,22);
for CountO = 1:x
    k1 = deltaNearFinal(CountO,1);
    k2 = deltaNearFinal(CountO,2);
    l1 = deltaNearFinal(CountO,3);
    l2 = deltaNearFinal(CountO,4);
    
    BalanceVect = BalanceVect .* 0;
    BalanceVect(1,k1) = BalanceVect(1,k1) -1;
    BalanceVect(1,k2) = BalanceVect(1,k2) -1;
    BalanceVect(1,l1) = BalanceVect(1,l1)+ 1;
    BalanceVect(1,l2) = BalanceVect(1,l2) +1;
    BalanceVect;
    InvestVect = InitialKet + BalanceVect;
    %%%% n values test
    
    if InvestVect(1,1) > 1
       % disp('mt');
    elseif any(InvestVect(1,12:22) > 1)
          % disp('n');
    else rr = rr +1;
        deltaFinal(rr,:) = deltaNearFinal(CountO,:);
    end
end



























