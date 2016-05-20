%%%%%% THIS ONLY WORKS ON BRAS!!!!!!
clc; clear all;
tic
A = [0 2 0 0 0 0 1 0 0 1 0 1 1 0 0 0 0 0 0 0];

KPos = [ 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20; %Matlab position
    0 0 1 0 1 0 1 0 1 0  1  0  1  0  1  0  1  0  1 0 ; %n value
    -1 0 0 1 1 2 2 3 3 4  4  5  5  6  6  7  7  8  8 9]; %Mt value

%A = [2 2];
z = 0;
firstsubtraction = zeros(1,length(A));
twosubtraction = zeros(1,length(A));
for a = 1:length(A)
    
    firstsubtraction = firstsubtraction.*0;
    firstsubtraction(1,a) = -1;
    B = A+firstsubtraction;
    for b = 1:length(A)
        
        twosubtraction = twosubtraction.*0;
        twosubtraction(1,b) = -1;
        C = B+twosubtraction;
        z = z + 1;
        MatA(z,:) = C;
        deltaP1(z,1) = a;
        deltaP1(z,2) = b;
    end
    
end

disp(['MAtA  is...', num2str(length(MatA))])
b = 0;
for CountB = 1:length(MatA)
    if any(MatA(CountB,:) < 0)
    else b = b + 1;
        MatB(b,:) = MatA(CountB,:);
        deltaP2(b,:) = deltaP1(CountB,:);
    end
end
disp(['MAtB  is...', num2str(length(MatB))])
MatB;
firstadd = zeros(1,length(A));
twoadd = zeros(1,length(A));
zzz=0;
[m n] = size(MatB);
for AAA = 1:m
    initialState = MatB(AAA,:);
    k1 = deltaP2(AAA,1);
    k2 = deltaP2(AAA,2);
    for BBB = 1:length(initialState)
        firstadd = firstadd.*0;
        firstadd(1,BBB) = 1;
        CCC = initialState + firstadd;
        for DDD = 1:length(initialState)
            twoadd = twoadd.*0;
            twoadd(1,DDD) = 1;
            EEE = CCC + twoadd;
            zzz = zzz + 1;
            MatC(zzz,:) = EEE;
            delta(zzz,1) = k1;
            delta(zzz,2) = k2;
            delta(zzz,3) = BBB;
            delta(zzz,4) = DDD;
        end
    end
    
end
disp(['MAtC  is...', num2str(length(MatC))])

L = 0;
for CountL = 1:length(delta)
    k1 = delta(CountL,1);
    k2 = delta(CountL,2);
    l1 = delta(CountL,3);
    l2 = delta(CountL,4);
    
    delta(CountL,5)  = KPos(2,k1);
    delta(CountL,6) = KPos(2,k2);
    delta(CountL,7) = KPos(2,l1);
    delta(CountL,8) = KPos(2,l2);
    
    delta(CountL,9)  = KPos(3,k1);
    delta(CountL,10) = KPos(3,k2);
    delta(CountL,11) = KPos(3,l1);
    delta(CountL,12) = KPos(3,l2); 
end
m = 0;
for CountM = 1:length(delta)
    
    if (delta(CountM,9) + delta(CountM,10)) == (delta(CountM,11) + delta(CountM,12))
        m = m+1;
        deltaMt(m,:) = delta(CountM,:);
    end
end

n=0;
[x y] = size(deltaMt);
for CountN = 1:x
    if ((deltaMt(CountN,5)+deltaMt(CountN,6)) < 2) && ((deltaMt(CountN,7)+deltaMt(CountN,8)) < 2)
        n = n + 1;
        deltaNMt(n,:) = deltaMt(CountN,:);
    end
end
x=0; y=0;
p = 0;
[x y] = size(deltaNMt);
for CountP = 1:x
    if ((deltaNMt(CountP,9) == -1) && (deltaNMt(CountP,10) == -1)) ||...
            ((deltaNMt(CountP,11) == -1) && (deltaNMt(CountP,12) == -1))
    else p = p + 1;
        deltaNMtOne(p,:) = deltaNMt(CountP,:);
    end
end
x=0; y=0;
q=0;
[x y] = size(deltaNMtOne);
for CountQ = 1:x
    if (((deltaNMtOne(CountQ,5) == 1) || (deltaNMtOne(CountQ,6) == 1)) && ((deltaNMtOne(CountQ,9) == -1) || (deltaNMtOne(CountQ,10) == -1)))...
            || (((deltaNMtOne(CountQ,7) == 1) || (deltaNMtOne(CountQ,8) == 1)) && ((deltaNMtOne(CountQ,11) == -1) || (deltaNMtOne(CountQ,12) == -1)))
    else q = q +1;
        deltaNMtOneNoMtLim(q,:) = deltaNMtOne(CountQ,:);
    end
    
end
x=0; y=0;
r = 0;
[x y] = size(deltaNMtOneNoMtLim);
for CountR = 1:x
    if (deltaNMtOneNoMtLim(CountR,11) + deltaNMtOneNoMtLim(CountR,12)) < 9
        r = r + 1;
        deltaFinal(r,:) = deltaNMtOneNoMtLim(CountR,:);
    end
end
deltaFinal


toc
