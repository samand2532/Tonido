clc; clear all;
tic
KPos = [ 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19; %Matlab position
          0 0 1 0 1 0 1 0 1 0  1  0  1  0  1  0  1  0  1 ; %n value
         -1 0 0 1 1 2 2 3 3 4  4  5  5  6  6  7  7  8  8]; %Mt value
     
     
    

Ca = 0;
for k1 = 1:18
    for k2 = 1:18
        for k3 = 1:18
            for k4 = 1:18
                Ca = Ca + 1;
                
                AllK(Ca,1) = k1;
                AllK(Ca,2) = k2;
                AllK(Ca,3) = k3;
                AllK(Ca,4) = k4;
            end
        end
    end
end
MatA = zeros(length(AllK),13);
Cb = 0;
for CountA = 1:length(AllK)
    Cb = Cb+1;
    k1 = AllK(Cb,1);
    k2 = AllK(Cb,2);
    k3 = AllK(Cb,3);
    k4 = AllK(Cb,4);
    
    mk1 = KPos(3,k1);
    mk2 = KPos(3,k2);
    ml1 = KPos(3,k3);
    ml2 = KPos(3,k4);
    
    nk1 = KPos(2,k1);
    nk2 = KPos(2,k2);
    nl1 = KPos(2,k3);
    nl2 = KPos(2,k4);
    
    MatA(Cb,1) = k1;
    MatA(Cb,2) = k2;
    MatA(Cb,3) = k3;
    MatA(Cb,4) = k4;
    
    MatA(Cb,5) = mk1;
    MatA(Cb,6) = mk2;
    MatA(Cb,7) = ml1;
    MatA(Cb,8) = ml2;
    
    MatA(Cb,9) = nk1;
    MatA(Cb,10) = nk2;
    MatA(Cb,11) = nl1;
    MatA(Cb,12) = nl2;
    
end

RemoveMt = (MatA(:,1:4) == 1);
SumRemoveMt = sum(RemoveMt,2);    
Cc = 0;
for CountB = 1:length(SumRemoveMt)
    
    if SumRemoveMt(CountB) < 2
    Cc = Cc + 1;
    MatB(Cc,1) = MatA(CountB,1);
    MatB(Cc,2) = MatA(CountB,2);
    MatB(Cc,3) = MatA(CountB,3);
    MatB(Cc,4) = MatA(CountB,4);
    
    MatB(Cc,5) = MatA(CountB,5);
    MatB(Cc,6) = MatA(CountB,6);
    MatB(Cc,7) = MatA(CountB,7);
    MatB(Cc,8) = MatA(CountB,8);
    
    MatB(Cc,9) = MatA(CountB,9);
    MatB(Cc,10) = MatA(CountB,10);
    MatB(Cc,11) = MatA(CountB,11);
    MatB(Cc,12) = MatA(CountB,12);
    end
end

RemoveN = (MatB(:,9:12) == 1);
SumRemoveN = sum(RemoveN,2);

Cd = 0;
for CountC = 1:length(SumRemoveN)
    
    if SumRemoveN(CountC) < 2
        Cd = Cd + 1;
        MatC(Cd,1) = MatB(CountC,1);
        MatC(Cd,2) = MatB(CountC,2);
        MatC(Cd,3) = MatB(CountC,3);
        MatC(Cd,4) = MatB(CountC,4);
        
        MatC(Cd,5) = MatB(CountC,5);
        MatC(Cd,6) = MatB(CountC,6);
        MatC(Cd,7) = MatB(CountC,7);
        MatC(Cd,8) = MatB(CountC,8);
        
        MatC(Cd,9) = MatB(CountC,9);
        MatC(Cd,10) = MatB(CountC,10);
        MatC(Cd,11) = MatB(CountC,11);
        MatC(Cd,12) = MatB(CountC,12);
    end
end

RemoveNa = MatC(:,9:12) == 1;
RemoveMt = MatC(:,1:4) == 1;

RemoveNaMt = [RemoveNa RemoveMt];
SumRemoveNaMt = sum(RemoveNaMt,2);

Ce = 0;
for CountD = 1:length(SumRemoveNaMt)
    if SumRemoveNaMt(CountD) < 2
        Ce = Ce + 1;  
        MatD(Ce,:) = MatC(CountD,:);       
    end
end

Cf = 0;
for CountE = 1:length(MatD)
    if (MatD(CountE,5) + MatD(CountE,6)) == (MatD(CountE,7) + MatD(CountE,8))
        Cf = Cf + 1;
        MatE(Cf,:) = MatD(CountE,:);
    end
end

Cg = 0;
for CountF = 1:length(MatE)
    if (MatE(CountF,5) + MatE(CountF,6)) <= 8
        Cg = Cg + 1;
        MatF(Cg,:) = MatE(CountF,:);
    end
end
    
for CountG = 1:length(MatF)
    CountG
    k1 = MatF(CountG,1);
    k2 = MatF(CountG,2);
    k3 = MatF(CountG,3);
    k4 = MatF(CountG,4);
    mk1 = MatF(CountG,5);
    mk2 = MatF(CountG,6);
    ml1 = MatF(CountG,7);
    ml2 = MatF(CountG,8);
    nk1 = MatF(CountG,9);
    nk2 = MatF(CountG,10);
    ml1 = MatF(CountG,11);
    ml2 = MatF(CountG,12);
    %%%%%%%% METHOD A %%%%%%%%%%%%%%
    Mt = abs(mk1)+abs(mk2)+abs(ml1)+abs(ml2);

    MomTerm = 1/(2^(Mt/2));
    
    fun = @(x) exp(-x).*(x.^(Mt/2)).*laguerreL(nk1,abs(mk1),(x/2)).*laguerreL(nk2,abs(mk2),(x/2))...
        .*laguerreL(nl1,abs(ml1),(x/2)).*laguerreL(nl2,abs(ml2),(x/2));
    
    I2func = integral(fun,0,inf);
    
    SqrtTerm = sqrt( 1/((factorial(nk1+abs(mk1)))*(factorial(nk2+abs(mk2)))*...
        (factorial(nl1+abs(ml1)))*(factorial(nl2+abs(ml2)))));
    
    Const = MomTerm*I2func*SqrtTerm;
    
    
    
    
    
    MatF(CountG,13) = Const;
    MatF(CountG,14) = MomTerm;
    MatF(CountG,15) = I2func;
    MatF(CountG,16) = SqrtTerm;
    
    
    
end
csvwrite('DeltaIntTry3LL2.csv',MatF);



toc
    
    
    
    
    
    
    
                