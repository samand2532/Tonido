%%%% new AnisoDelta that include n=1 to n=1 transition
clc; clear all;

tic

KPos = [ 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19; %Matlab position
    0 0 1 0 1 0 1 0 1 0  1  0  1  0  1  0  1  0  1 ; %n value
    -1 0 0 1 1 2 2 3 3 4  4  5  5  6  6  7  7  8  8]; %Mt value

ConstA = 0;
for k1 = 1:19
    for k2 = 1:19
        
        ConstA = ConstA + 1;
        MatA(ConstA,1) = k1;
        MatA(ConstA,2) = k2;
        
        MatA(ConstA,3) = KPos(2,k1);
        MatA(ConstA,4) = KPos(2,k2);
        
        MatA(ConstA,5) = KPos(3,k1);
        MatA(ConstA,6) = KPos(3,k2);
        
        
    end
end

ConstB = 0;
for CountB = 1:length(MatA)
    if MatA(CountB,6) == (MatA(CountB,5) + 2)
        ConstB = ConstB + 1;
        MatB(ConstB,:) = MatA(CountB,:);
    elseif MatA(CountB,6) == (MatA(CountB,5) - 2)
        ConstB = ConstB + 1;
        MatB(ConstB,:) = MatA(CountB,:);
    end
end


MatB
toc
