%%% Laguerre Values

clc;clear all;
tic
LengthVect =19;
MomValue = [ -1 0 0 1 1 2 2 3 3 4 4 5 5 6 6 7 7 8 8];  

KPos = [ 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19; %Matlab position
         0 0 1 0 1 0 1 0 1 0  1  0  1  0  1  0  1  0  1 ; %n value
        -1 0 0 1 1 2 2 3 3 4  4  5  5  6  6  7  7  8  8]; %Mt value

IntDelta = csvread('IntDelta.csv');
PotDelta = csvread('PotDelta.csv');



%%%%% nk1 | nk2 | nl1 | nl2| mk1| mk2| ml1| ml2| Lag


%%%%% Int Term
LagTable = zeros(length(IntDelta),9);
for IntCount = 1:length(IntDelta)
    IntCount
    k1 = IntDelta(IntCount,1);
    k2 = IntDelta(IntCount,2);
    l1 = IntDelta(IntCount,3);
    l2 = IntDelta(IntCount,4);
    
    mk1 = KPos(3,k1);
    mk2 = KPos(3,k2);
    ml1 = KPos(3,l1);
    ml2 = KPos(3,l2);
    nk1 = KPos(2,k1);
    nk2 = KPos(2,k2);
    nl1 = KPos(2,l1);
    nl2 = KPos(2,l2);
    
    LagTable(IntCount,1) = nk1;
    LagTable(IntCount,2) = nk2;
    LagTable(IntCount,3) = nl1;
    LagTable(IntCount,4) = nl2;
    LagTable(IntCount,5) = mk1;
    LagTable(IntCount,6) = mk2;
    LagTable(IntCount,7) = ml1;
    LagTable(IntCount,8) = ml2;
    
     Mt = (abs(mk1)+abs(mk2)+abs(ml1)+abs(ml2))/2;
     fun = @(x) exp(-x).*(x.^(Mt)).* (laguerreL(nk1,abs(mk1),(x/2))) .* (laguerreL(nk2,abs(mk2),(x/2))) .* (laguerreL(nl1,abs(ml1),(x/2))) .* (laguerreL(nl2,abs(ml2),(x/2)));
     I2func = integral(fun,0,inf);
     
     LagTable(IntCount,9) = I2func;
end
toc