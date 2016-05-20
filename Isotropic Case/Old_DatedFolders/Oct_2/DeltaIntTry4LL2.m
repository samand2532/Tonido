clc; clear all;
tic

KPos = [ 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19; %Matlab position
          0 0 1 0 1 0 1 0 1 0  1  0  1  0  1  0  1  0  1 ; %n value
         -1 0 0 1 1 2 2 3 3 4  4  5  5  6  6  7  7  8  8]; %Mt value
     
     Ca = 0;
     
     for k1 = 1:18;
         for k2 = 1:18;
             for k3 = 1:18;
                 for k4 = 1:18;
                     Ca = Ca + 1;
                     
                     AllK(Ca,1) = k1;
                     AllK(Ca,2) = k2;
                     AllK(Ca,3) = k3;
                     AllK(Ca,4) = k4;
                     
                 end
             end
         end
     end
     MatA = zeros(length(AllK),17);
     for CountA = 1:length(AllK)
         
         k1 = AllK(CountA,1);
         k2 = AllK(CountA,2);
         k3 = AllK(CountA,3);
         k4 = AllK(CountA,4);
         
         mk1 = KPos(3,k1);
         mk2 = KPos(3,k2);
         ml1 = KPos(3,k3);
         ml2 = KPos(3,k4);
         
         nk1 = KPos(2,k1);
         nk2 = KPos(2,k2);
         nl1 = KPos(2,k3);
         nl2 = KPos(2,k4);
         
         MatA(CountA,1) = k1;
         MatA(CountA,2) = k2;
         MatA(CountA,3) = k3;
         MatA(CountA,4) = k4;
         
         MatA(CountA,5) = mk1;
         MatA(CountA,6) = mk2;
         MatA(CountA,7) = ml1;
         MatA(CountA,8) = ml2;
         
         MatA(CountA,9) = nk1;
         MatA(CountA,10) = nk2;
         MatA(CountA,11) = nl1;
         MatA(CountA,12) = nl2;
         
     end
         
         RemoveMt = (MatA(:,5:8) == -1);
         sumRemoveMt = sum(RemoveMt,2);
         
         Cb = 0;
         for CountB = 1:length(MatA)
             if sumRemoveMt(CountB) < 2
                 Cb = Cb + 1;
                 MatB(Cb,:) = MatA(CountB,:);
             end
         end
       toc  
       
       Cc = 0;
       for CountC = 1:length(MatB)
           mL = MatB(CountC,5) + MatB(CountC,6);
           mR = MatB(CountC,7) + MatB(CountC,8);
           if mL == mR
               Cc = Cc + 1;
               MatC(Cc,:) = MatB(CountC,:);
           end
       end
       
       Cd = 0;
       for CountD = 1:length(MatC)
           mL = MatC(CountD,5) + MatC(CountD,6);
           if mL < 9
               Cd = Cd + 1;
               MatD(Cd,:) = MatC(CountD,:);
           end
       end
       
       RemoveN = (MatD(:,9:12) == 1);
       sumRemoveN = sum(RemoveN,2);
       
       Ce = 0;
       for CountE = 1:length(MatD)
           if sumRemoveN(CountE) < 2
               Ce = Ce + 1;
               MatE(Ce,:) = MatD(CountE,:);
           end
       end
       
       %%%%% SEction removes occurance of Mt = -1 and another n=1
       
       Mtlogic = (MatE(:,5:8) == -1);
       Nlogic = (MatE(:,9:12) == 1);
       MtNlogic = [Mtlogic Nlogic];
       sumMtNlogic = sum(MtNlogic,2);
       
       Cf = 0;
       for CountF = 1:length(MatE)
           if sumMtNlogic(CountF) < 2
               Cf = Cf + 1;
               MatF(Cf,:) = MatE(CountF,:);
           end
       end
       
         %%% Of Form
         %%% k1|k2|k3|k4|mk1|mk2|ml1|ml2|nk1|nk2|nl1|nl2|Momterm|RootTerm|I2|Const|ConstforLLL
         %%% 
         %%% CONSTforLLL is only valid when n_1:4 == 0
         %%%  1|2 | 3| 4| 5 | 6 | 7 | 8 | 9 | 10|11 |12 | 13    |14      |15|   16  |  | 
         
         
        for CountG = 1:length(MatF)
            CountG
            mk1 = MatF(CountG,5);
            mk2 = MatF(CountG,6);
            ml1 = MatF(CountG,7);
            ml2 = MatF(CountG,8);
            
            nk1 = MatF(CountG,9);
            nk2 = MatF(CountG,10);
            nl1 = MatF(CountG,11);
            nl2 = MatF(CountG,12);
            
            Mt = abs(mk1)+abs(mk2)+abs(ml1)+abs(ml2);
            MomTerm = 1/(2^((Mt)/2));
            MatF(CountG,13) = MomTerm;
            
            SqrtTermTop = factorial(nk1)*factorial(nk2)*factorial(nl1)*factorial(nl2);
            SqrtTermBot = factorial(nk1+abs(mk1))*factorial(nk2+abs(mk2))*factorial(nl1+abs(ml1))*factorial(nl2+abs(ml2));
            
            SqrtTerm = SqrtTermTop/SqrtTermBot; 
            MatF(CountG,14) = SqrtTerm;
            
            amk1 = abs(mk1); amk2 = abs(mk2); aml1 = abs(ml1); aml2 = abs(ml2);
            fun = @(x) exp(-x).*x.^(Mt/2).*laguerreL(nk1,amk1,(x/2)).*laguerreL(nk2,amk2,(x/2)).*laguerreL(nl1,aml1,(x/2)).*laguerreL(nl2,aml2,(x/2));
            I2func = integral(fun,0,inf);
            MatF(CountG,15) = I2func;
            
            I2a = factorial((amk1 + amk2 + aml1 + aml2)/2);
            
            Comb = MomTerm.*SqrtTerm.*I2func;
            Comba = MomTerm.*SqrtTerm.*I2a;                                                                  
            MatF(CountG,16) = Comb;
            MatF(CountG,17) = Comba;
        end
        
        csvwrite('DeltaIntLL2Try4.csv',MatF);
        
        
         
         
         toc
         
         
         
     
                     