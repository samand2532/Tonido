clc; clear all;
tic
KPos = [ 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19; %Matlab position
         0 0 1 0 1 0 1 0 1 0  1  0  1  0  1  0  1  0  1 ; %n value
        -1 0 0 1 1 2 2 3 3 4  4  5  5  6  6  7  7  8  8]; %Mt value
LLVal = 2; %1 = LLL 2 = LL2


if LLVal == 1
    %%%%%%%% WRONG! BALANCED K RATHER THAN MT!!!!
    Ua = 0; Ub = 0; Uc = 0; Ud = 0;
    for d1 = 2:2:18
        for d2 = 2:2:18
            for d3 = 2:2:18
                for d4 = 2:2:18
                    Ua = Ua + 1;
                    Mat_A(Ua,1) = d1;
                    Mat_A(Ua,2) = d2;
                    Mat_A(Ua,3) = d3;
                    Mat_A(Ua,4) = d4;
                    
                    if (Mat_A(Ua,1)+Mat_A(Ua,2)) == (Mat_A(Ua,3)+Mat_A(Ua,4))
                        Ub = Ub + 1;
                        Mat_B(Ub,:) = Mat_A(Ua,:);
                    end
                end
            end
        end
    end
    
    
    for CountA = 1:length(Mat_B)
        k1 = Mat_B(CountA,1);
        k2 = Mat_B(CountA,2);
        Mtk1 = KPos(3,k1);
        Mtk2 = KPos(3,k2);
        if (Mtk1+Mtk2) <= 8
            Uc = Uc +1;
            Mat_C(Uc,:) = Mat_B(CountA,:);
        end
    end
    
    for CountB = 1:length(Mat_C)
        Ud = Ud + 1;
        
        k1 = Mat_C(Ud,1);
        k2 = Mat_C(Ud,2);
        l1 = Mat_C(Ud,3);
        l2 = Mat_C(Ud,4);
        
        mk1 = KPos(3,k1);
        mk2 = KPos(3,k2);
        ml1 = KPos(3,l1);
        ml2 = KPos(3,l2);
        
        nk1 = KPos(2,k1);
        nk2 = KPos(2,k2);
        nl1 = KPos(2,l1);
        nl2 = KPos(2,l2);
        
        FLT(Ud,1) = nk1;
        FLT(Ud,2) = nk2;
        FLT(Ud,3) = nl1;
        FLT(Ud,4) = nl2;
        FLT(Ud,5) = mk1;
        FLT(Ud,6) = mk2;
        FLT(Ud,7) = ml1;
        FLT(Ud,8) = ml2;
        
        Mt = (abs(mk1)+abs(mk2)+abs(ml1)+abs(ml2));
        MomTerm = 1/(2^(Mt/2));
        
        fun = @(x) exp(-x).*(x.^(Mt/2)).* (laguerreL(nk1,abs(mk1)).*(x/2)) .* (laguerreL(nk2,abs(mk2)).*(x/2)) .* (laguerreL(nl1,abs(ml1)).*(x/2)) .* (laguerreL(nl2,abs(ml2)).*(x/2));
        I2func = integral(fun,0,inf);
        
        PiTerm = factorial(nk1)*factorial(nk2)*factorial(nl1)*factorial(nl2)/...
            factorial(nk1+abs(mk1))*factorial(nk2+abs(mk2))*factorial(nl1+abs(ml1))*factorial(nl2+abs(ml2));
        
        FLT(Ud,9) = I2func*MomTerm*PiTerm;
        csvwrite('FLT_LLL.csv',FLT);
        
    end
end

if LLVal == 2
    
    Ua = 0; Ub = 0; Uc = 0; Ud = 0;  Ue = 0; Uf = 0;
    for d1 = 1:18
        for d2 = 1:18
            for d3 = 1:18
                for d4 = 1:18
                    Ua = Ua + 1;
                    Mat_A(Ua,1) = d1;
                    Mat_A(Ua,2) = d2;
                    Mat_A(Ua,3) = d3;
                    Mat_A(Ua,4) = d4;
                    
                    if (Mat_A(Ua,1)+Mat_A(Ua,2)) == (Mat_A(Ua,3)+Mat_A(Ua,4))
                        Ub = Ub + 1;
                        Mat_B(Ub,:) = Mat_A(Ua,:);
                    end
                end
            end
        end
    end

    for CountA = 1:length(Mat_B)
        k1 = Mat_B(CountA,1);
        k2 = Mat_B(CountA,2);
        Mtk1 = KPos(3,k1);
        Mtk2 = KPos(3,k2);
        if (Mtk1+Mtk2) <= 8
            Uc = Uc +1;
            Mat_C(Uc,:) = Mat_B(CountA,:);
        end
    end
    
   %%% removes multiple n's 
   for CountB = 1:length(Mat_C)
       k1 = Mat_C(CountB,1);
       k2 = Mat_C(CountB,2);
       l1 = Mat_C(CountB,3);
       l2 = Mat_C(CountB,4);
       % TO DO PULL OUT N VALUES AND DESELECT MORE THAN 
       %1 occurance, then do the same for multiple -1MT.
       Mat_D(CountB,1) = KPos(2,k1);
       Mat_D(CountB,2) = KPos(2,k2);
       Mat_D(CountB,3) = KPos(2,l1);
       Mat_D(CountB,4) = KPos(2,l2);
       
       
       
    
    
   end
    Array_A = (Mat_D == 1);
    Array_B = sum(Array_A,2);
    
    for CountC = 1:length(Array_B)
        if Array_B(CountC)  <= 1
            Ud=Ud+1;
            Mat_E(Ud,:) = Mat_C(CountC,:);
        end
    end
    
    %%%% remove multiple Mt = -1 values
    
    Array_C = (Mat_E == 1);
    Array_D = sum(Array_C,2);
    
    for CountD = 1:length(Array_D)
        if Array_D(CountD) < 2
            Ue = Ue + 1;
            Mat_F(Ue,:) = Mat_E(CountD,:);
        end
    end
    
    for CountE = 1:length(Mat_F)
        Uf = Uf + 1;
        
        k1 = Mat_F(Uf,1);
        k2 = Mat_F(Uf,2);
        l1 = Mat_F(Uf,3);
        l2 = Mat_F(Uf,4);
        
        mk1 = KPos(3,k1);
        mk2 = KPos(3,k2);
        ml1 = KPos(3,l1);
        ml2 = KPos(3,l2);
        
        nk1 = KPos(2,k1);
        nk2 = KPos(2,k2);
        nl1 = KPos(2,l1);
        nl2 = KPos(2,l2);
        
        FLT(Uf,1) = nk1;
        FLT(Uf,2) = nk2;
        FLT(Uf,3) = nl1;
        FLT(Uf,4) = nl2;
        FLT(Uf,5) = mk1;
        FLT(Uf,6) = mk2;
        FLT(Uf,7) = ml1;
        FLT(Uf,8) = ml2;
        
        Mt = (abs(mk1)+abs(mk2)+abs(ml1)+abs(ml2));
        MomTerm = 1/(2^(Mt/2));
        
        fun = @(x) exp(-x).*(x.^(Mt/2)).* (laguerreL(nk1,abs(mk1)).*(x/2)) .* (laguerreL(nk2,abs(mk2)).*(x/2)) .* (laguerreL(nl1,abs(ml1)).*(x/2)) .* (laguerreL(nl2,abs(ml2)).*(x/2));
        I2func = integral(fun,0,inf);
        
        PiTerm = factorial(nk1)*factorial(nk2)*factorial(nl1)*factorial(nl2)/...
            factorial(nk1+abs(mk1))*factorial(nk2+abs(mk2))*factorial(nl1+abs(ml1))*factorial(nl2+abs(ml2));
        
        FLT(Ud,9) = I2func*MomTerm*PiTerm;
    end
    
    csvwrite('FLT_LL2.csv',FLT);
end
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
                        
                    