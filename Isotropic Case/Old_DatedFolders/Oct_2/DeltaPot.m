clc; clear all;

KPos = [ 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19; %Matlab position
    0 0 1 0 1 0 1 0 1 0  1  0  1  0  1  0  1  0  1 ; %n value
    -1 0 0 1 1 2 2 3 3 4  4  5  5  6  6  7  7  8  8]; %Mt value


g = 1; N = 6; A = 0.00; LLVal = 1; % 1 = LLL, 2 = LL2

if LLVal == 1;
   Ca = 0;
    for k1 = 2:2:18;
        for k2 = 2:2:18;
            Ca = Ca + 1;
            MatA(Ca,1) = k1;
            MatA(Ca,2) = k2;
            
        end
    end
    
    Cb = 0;
    for CountA = 1:length(MatA)
        Cb = Cb+1;
        
        nk1 = MatA(CountA,1);
        nk2 = MatA(CountA,2);
        
        NMat(Cb,1) = KPos(2,nk1);
        NMat(Cb,2) = KPos(2,nk2);
    end
    
    NPos = (NMat == 1);
    NSum = sum(NPos,2);
    
    Cc  =0;
    for CountB = 1:length(NSum)
        
        if NSum(CountB) < 2
            Cc = Cc + 1;
            
            MatB(Cc,1) = MatA(CountB,1);
            MatB(Cc,2) = MatA(CountB,2);
            
            
        end
    end
    
    %%%% Remove multipe Mts = -1
    
    RemoveMt1 = (MatB == 1);
    SumRemoveMt1 = sum(RemoveMt1,2);
    Cd = 0;
    for CountC = 1:length(SumRemoveMt1)
        
        if SumRemoveMt1(CountC) < 2
            Cd = Cd + 1;
            MatC(Cd,1) = MatB(CountC,1);
            MatC(Cd,2) = MatB(CountC,2);
            
        end
    end
    
    Ce = 0;Cf = 0;
    for CountD = 1:length(MatC)
        k1 = MatC(CountD,1);
        k2 = MatC(CountD,2);
        
        mk1 = KPos(3,k1);
        mk2 = KPos(3,k2);
        
        nk1 = KPos(2,k1);
        nk2 = KPos(2,k2);
        
        if mk1 == (mk2 + 2);
            Ce = Ce + 1;
            
            MatDPlus(Ce,1) = MatC(CountD,1);
            MatDPlus(Ce,2) = MatC(CountD,2);
            MatDPlus(Ce,3) = nk1;
            MatDPlus(Ce,4) = nk2;
            MatDPlus(Ce,5) = mk1;
            MatDPlus(Ce,6) = mk2;
            
            
            
        elseif mk1 == (mk2 - 2);
            Cf = Cf +1;
            
            MatDMinus(Cf,1) = MatC(CountD,1);
            MatDMinus(Cf,2) = MatC(CountD,2);
            MatDMinus(Cf,3) = nk1;
            MatDMinus(Cf,4) = nk2;
            MatDMinus(Cf,5) = mk1;
            MatDMinus(Cf,6) = mk2;
            
        end
    end
    
    delta = [MatDMinus; MatDPlus];
    
    for CountE = 1:length(delta)
        
        k1 = delta(CountE,1);
        k2 = delta(CountE,2);
        nk1 = delta(CountE,3);
        nk2 = delta(CountE,4);
        mk1 = delta(CountE,5);
        mk2 = delta(CountE,6);
        
        RootTerm = sqrt((factorial(nk1)*factorial(nk2))/((factorial(nk1 + abs(mk1)))*factorial(nk2 + abs(mk2))));
        
        Mt = (abs(mk1)+abs(mk2)+2)/2;
        fun = @(x) exp(-x).*(x.^(Mt)).* (laguerreL(nk1,abs(mk1),(x))) .* (laguerreL(nk2,abs(mk2),(x)));
        I1func = integral(fun,0,inf);
        
        TotalConst = RootTerm*I1func;
       
        
        LagTable(CountE,1) = k1;
        LagTable(CountE,2) = k2;
        LagTable(CountE,3) = nk1;
        LagTable(CountE,4) = nk2;
        LagTable(CountE,5) = mk1;
        LagTable(CountE,6) = mk2;
        LagTable(CountE,7) = TotalConst;
        
    end
    
    csvwrite('LagPotLLL.csv',LagTable)
    
    
end



if LLVal == 2;
    Ca = 0;
    for k1 = 1:18;
        for k2 = 1:18;
            Ca = Ca + 1;
            MatA(Ca,1) = k1;
            MatA(Ca,2) = k2;
            
        end
    end
    
    %%%% remove multiple n = 1's
    Cb = 0;
    for CountA = 1:length(MatA)
        Cb = Cb+1;
        
        nk1 = MatA(CountA,1);
        nk2 = MatA(CountA,2);
        
        NMat(Cb,1) = KPos(2,nk1);
        NMat(Cb,2) = KPos(2,nk2);
    end
    
    NPos = (NMat == 1);
    NSum = sum(NPos,2);
    
    Cc  =0;
    for CountB = 1:length(NSum)
        
        if NSum(CountB) < 2
            Cc = Cc + 1;
            
            MatB(Cc,1) = MatA(CountB,1);
            MatB(Cc,2) = MatA(CountB,2);
            
            
        end
    end
    
    %%%% Remove multipe Mts = -1
    
    RemoveMt1 = (MatB == 1);
    SumRemoveMt1 = sum(RemoveMt1,2);
    Cd = 0;
    for CountC = 1:length(SumRemoveMt1)
        
        if SumRemoveMt1(CountC) < 2
            Cd = Cd + 1;
            MatC(Cd,1) = MatB(CountC,1);
            MatC(Cd,2) = MatB(CountC,2);
            
        end
    end
    
    Ce = 0;Cf = 0;
    for CountD = 1:length(MatC)
        k1 = MatC(CountD,1);
        k2 = MatC(CountD,2);
        
        mk1 = KPos(3,k1);
        mk2 = KPos(3,k2);
        
        nk1 = KPos(2,k1);
        nk2 = KPos(2,k2);
        
        if mk1 == (mk2 + 2);
            Ce = Ce + 1;
            
            MatDPlus(Ce,1) = MatC(CountD,1);
            MatDPlus(Ce,2) = MatC(CountD,2);
            MatDPlus(Ce,3) = nk1;
            MatDPlus(Ce,4) = nk2;
            MatDPlus(Ce,5) = mk1;
            MatDPlus(Ce,6) = mk2;
            
            
            
        elseif mk1 == (mk2 - 2);
            Cf = Cf +1;
            
            MatDMinus(Cf,1) = MatC(CountD,1);
            MatDMinus(Cf,2) = MatC(CountD,2);
            MatDMinus(Cf,3) = nk1;
            MatDMinus(Cf,4) = nk2;
            MatDMinus(Cf,5) = mk1;
            MatDMinus(Cf,6) = mk2;
            
        end
    end
    
    delta = [MatDMinus; MatDPlus];
    
    for CountE = 1:length(delta)
        
        k1 = delta(CountE,1);
        k2 = delta(CountE,2);
        nk1 = delta(CountE,3);
        nk2 = delta(CountE,4);
        mk1 = delta(CountE,5);
        mk2 = delta(CountE,6);
        
        RootTerm = sqrt((factorial(nk1)*factorial(nk2))/((factorial(nk1 + abs(mk1)))*factorial(nk2 + abs(mk2))));
        
        Mt = (abs(mk1)+abs(mk2)+2)/2;
        fun = @(x) exp(-x).*(x.^(Mt)).* (laguerreL(nk1,abs(mk1),(x))) .* (laguerreL(nk2,abs(mk2),(x)));
        I1func = integral(fun,0,inf);
        
        TotalConst = RootTerm*I1func;
       
        
        LagTable(CountE,1) = k1;
        LagTable(CountE,2) = k2;
        LagTable(CountE,3) = nk1;
        LagTable(CountE,4) = nk2;
        LagTable(CountE,5) = mk1;
        LagTable(CountE,6) = mk2;
        LagTable(CountE,7) = TotalConst;
        
    end
        csvwrite('LagPotLL2.csv',LagTable)
    end
    