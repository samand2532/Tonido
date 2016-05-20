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
%%%% Remove multiple n's
nRem = MatA(:,3:4) == 1;
sumnRem = sum(nRem,2);

ConstB = 0;
for CountB = 1:length(sumnRem)
    if sumnRem(CountB) < 2
        ConstB = ConstB + 1;
        MatB(ConstB,:) = MatA(CountB,:);
    end
end

%%%% Remove multiple Mt = -1's

mRem = MatB(:,5:6) == -1;
summRem = sum(mRem,2);

ConstC = 0;
for CountC = 1:length(MatB)
    if summRem(CountC) < 2
        ConstC = ConstC + 1;
        MatC(ConstC,:) = MatB(CountC,:);
    end
end

%%%%% Stop n and mt=-1 occuring simul...
n1rem = MatC(:,3:4) == 1;
m1rem = MatC(:,5:6) == -1;
rem1mat = [n1rem m1rem];
sumrem1mat = sum(rem1mat,2);

ConstD = 0;
for CountD = 1:length(sumrem1mat)
    if sumrem1mat(CountD) < 2
        ConstD = ConstD + 1;
        MatD(ConstD,:) = MatC(CountD,:);
    end
end

%%%%%%% mk2 +/- mk1

ConstE = 0;
for CountE = 1:length(MatD)
    if MatD(CountE,5) == (MatD(CountE,6) + 2)
        ConstE = ConstE + 1;
        MatE(ConstE,:) = MatD(CountE,:);
    elseif MatD(CountE,5) == (MatD(CountE,6) - 2)
        ConstE = ConstE + 1;
        MatE(ConstE,:) = MatD(CountE,:);
    end
end
%%% Limits to Mt <9
ConstF = 0;
for CountF = 1:length(MatE)
    if (MatE(CountF,5) + MatE(CountF,6)) < 9
        ConstF = ConstF + 1;
        MatF(ConstF,:) = MatE(CountF,:);
    end
end
MatF
toc

