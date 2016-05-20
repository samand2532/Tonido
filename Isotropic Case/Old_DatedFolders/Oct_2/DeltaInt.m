clc; clear all;

KPos = [ 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19; %Matlab position
          0 0 1 0 1 0 1 0 1 0  1  0  1  0  1  0  1  0  1 ; %n value
         -1 0 0 1 1 2 2 3 3 4  4  5  5  6  6  7  7  8  8]; %Mt value
     
     
     g = 1; N = 6; A = 0.00; LLVal = 1; % 1 = LLL, 2 = LL2
     
     
     if LLVal == 1;
         Ca = 0; 
         for k1 = 2:2:18
             for k2 = 2:2:18
                 for k3 = 2:2:18
                     for k4 = 2:2:18
                         Ca = Ca + 1;
                         MatA(Ca,1) = k1;
                         MatA(Ca,2) = k2;
                         MatA(Ca,3) = k3;
                         MatA(Ca,4) = k4;
                         
                     end
                 end
             end
         end
         Cb = 0;
         for CountA = 1:length(MatA)
             Cb = Cb + 1;
             mk1 = MatA(Cb,1);
             mk2 = MatA(Cb,2);
             ml1 = MatA(Cb,3);
             ml2 = MatA(Cb,4);
             
             MatB(Cb,1) = MatA(Cb,1);
             MatB(Cb,2) = MatA(Cb,2);
             MatB(Cb,3) = MatA(Cb,3);
             MatB(Cb,4) = MatA(Cb,4);
             MatB(Cb,5) = KPos(3,mk1);
             MatB(Cb,6) = KPos(3,mk2);
             MatB(Cb,7) = KPos(3,ml1);
             MatB(Cb,8) = KPos(3,ml2);
         end
         %%%%%%  section balances mk1+mk2 == ml1 + ml2
         Ce = 0;
         for CountE = 1:length(MatB)
             
             if (MatB(CountE,5)+MatB(CountE,6)) == (MatB(CountE,7)+MatB(CountE,8))
            
             Ce = Ce +1;
              MatC(Ce,1) = MatB(CountE,1);
               MatC(Ce,2) = MatB(CountE,2);
                MatC(Ce,3) = MatB(CountE,3);
                 MatC(Ce,4) = MatB(CountE,4);
             MatC(Ce,5) = MatB(CountE,5);
             MatC(Ce,6) = MatB(CountE,6);
             MatC(Ce,7) = MatB(CountE,7);
             MatC(Ce,8) = MatB(CountE,8);
             end
         end
         Cd = 0;
         
         for CountC = 1:length(MatC)
             if (MatC(CountC,5) + MatC(CountC,6)) <= 8
                 Cd = Cd +1;
                  MatD(Cd,1) = MatC(CountC,1);
                   MatD(Cd,2) = MatC(CountC,2);
                    MatD(Cd,3) = MatC(CountC,3);
                     MatD(Cd,4) = MatC(CountC,4);
                 MatD(Cd,5) = MatC(CountC,5);
                 MatD(Cd,6) = MatC(CountC,6);
                 MatD(Cd,7) = MatC(CountC,7);
                 MatD(Cd,8) = MatC(CountC,8);
             end
         end
         
         for CountD = 1:length(MatD)
         mk1 = MatD(CountD,5);
         mk2 = MatD(CountD,6);
         ml1 = MatD(CountD,7);
         ml2 = MatD(CountD,8);
        
         nk1 = 0;
         nk2 = 0;
         nl1 = 0;
         nl2 = 0;
         
         LagTable(CountD,1) = MatD(CountD,1);
         LagTable(CountD,2) = MatD(CountD,2);
         LagTable(CountD,3) = MatD(CountD,3);
         LagTable(CountD,4) = MatD(CountD,4);
         LagTable(CountD,5) = nk1;
         LagTable(CountD,6) = nk2;
         LagTable(CountD,7) = nl1;
         LagTable(CountD,8) = nl2;
         LagTable(CountD,9) = mk1;
         LagTable(CountD,10) = mk2;
         LagTable(CountD,11) = ml1;
         LagTable(CountD,12) = ml2;
         
         
         PiTerm = sqrt(1/ (factorial(abs(mk1))*factorial(abs(mk2))*factorial(abs(ml1))*factorial(abs(ml2))));
         
         Mt = (abs(mk1)+abs(mk2)+abs(ml1)+abs(ml2))/2;
         fun = @(x) exp(-x).*(x.^(Mt)).* (laguerreL(nk1,abs(mk1),(x/2))) .* (laguerreL(nk2,abs(mk2),(x/2))) .* (laguerreL(nl1,abs(ml1),(x/2))) .* (laguerreL(nl2,abs(ml2),(x/2)));
         
         I2func = integral(fun,0,inf);
         
         MomTerm = 1/(2^((abs(mk1)+abs(mk2)+abs(ml1)+abs(ml2))/2));
         
         TotalConst = PiTerm.*I2func.*MomTerm;
        
         LagTable(CountD,13) = TotalConst;
         end
   
         
         csvwrite('LagIntLLL.csv',LagTable)        
     end
     
     
     if LLVal == 2;
         Ca = 0; 
         for k1 = 1:18
             for k2 = 1:18
                 for k3 = 1:18
                     for k4 = 1:18
                         Ca = Ca + 1;
                         MatA(Ca,1) = k1;
                         MatA(Ca,2) = k2;
                         MatA(Ca,3) = k3;
                         MatA(Ca,4) = k4;
                         
                     end
                 end
             end
         end
         
         %%%% REmove mutliple N's
         Cb = 0;
         for CountA = 1:length(MatA)
             Cb = Cb + 1;
            
             nk1 = MatA(CountA,1);
             nk2 = MatA(CountA,2);
             nl1 = MatA(CountA,3);
             nl2 = MatA(CountA,4);
             
             NMat(Cb,1) = KPos(2,nk1);
             NMat(Cb,2) = KPos(2,nk2);
             NMat(Cb,3) = KPos(2,nl1);
             NMat(Cb,4) = KPos(2,nl2);
         end
         
         NPos = (NMat == 1);
         NSum = sum(NPos,2);
         Cc = 0;
         for CountB = 1:length(NSum)
             
             
             if NSum(CountB) <= 1
                 Cc = Cc + 1;
                 MatB(Cc,1) = MatA(CountB,1);
                 MatB(Cc,2) = MatA(CountB,2);
                 MatB(Cc,3) = MatA(CountB,3);
                 MatB(Cc,4) = MatA(CountB,4);
                 MatB(Cc,5) = NMat(CountB,1);
                 MatB(Cc,6) = NMat(CountB,2);
                 MatB(Cc,7) = NMat(CountB,3);
                 MatB(Cc,8) = NMat(CountB,4);
        
             end
         end
         
         
         %%%% Remove multipe Mts = -1
         
         RemoveMt1 = (MatB(:,1:4) == 1);
         SumRemoveMt1 = sum(RemoveMt1,2);
         Cd = 0;
         for CountC = 1:length(SumRemoveMt1)
             
             if SumRemoveMt1(CountC) < 2
             Cd = Cd + 1;
             MatC(Cd,1) = MatB(CountC,1);
             MatC(Cd,2) = MatB(CountC,2);
             MatC(Cd,3) = MatB(CountC,3);
             MatC(Cd,4) = MatB(CountC,4);
             MatC(Cd,5) = MatB(CountC,5);
             MatC(Cd,6) = MatB(CountC,6);
             MatC(Cd,7) = MatB(CountC,7);
             MatC(Cd,8) = MatB(CountC,8);
       
             end
         end
  
         %%%% Balance left and right.
         Ce = 0;
         for CountD = 1:length(MatC)
             k1 = MatC(CountD,1);
             k2 = MatC(CountD,2);
             l1 = MatC(CountD,3);
             l2 = MatC(CountD,4);
             
             mk1 = KPos(3,k1);
             mk2 = KPos(3,k2);
             ml1 = KPos(3,l1);
             ml2 = KPos(3,l2);
             
             if (mk1+mk2) == (ml1+ml2)
                 Ce = Ce + 1;
             
             
             MatD(Ce,1) = MatC(CountD,1);
             MatD(Ce,2) = MatC(CountD,2);
             MatD(Ce,3) = MatC(CountD,3);
             MatD(Ce,4) = MatC(CountD,4);
             MatD(Ce,5) = MatC(CountD,5);
             MatD(Ce,6) = MatC(CountD,6);
             MatD(Ce,7) = MatC(CountD,7);
             MatD(Ce,8) = MatC(CountD,8);
             MatD(Ce,9) = mk1;
             MatD(Ce,10) = mk2;
             MatD(Ce,11) = ml1;
             MatD(Ce,12) = ml2;
                          
             end
         end
     
         %%% Limit Mt <= 8 (per side)
         Cf = 0;
         for CountE  =1:length(MatD)
             k1 = MatD(CountE,1);
             k2 = MatD(CountE,2);
             mk1 = KPos(3,k1);
             mk2 = KPos(3,k2);
             
             if (mk1 + mk2) <= 8
                 Cf = Cf +1;
                 MatE(Cf,1) = MatD(CountE,1);
                 MatE(Cf,2) = MatD(CountE,2);
                 MatE(Cf,3) = MatD(CountE,3);
                 MatE(Cf,4) = MatD(CountE,4);
                 MatE(Cf,5) = MatD(CountE,5);
                 MatE(Cf,6) = MatD(CountE,6);
                 MatE(Cf,7) = MatD(CountE,7);
                 MatE(Cf,8) = MatD(CountE,8);
                 MatE(Cf,9) = MatD(CountE,9);
                 MatE(Cf,10) = MatD(CountE,10);
                 MatE(Cf,11) = MatD(CountE,11);
                 MatE(Cf,12) = MatD(CountE,12);
             end
         end
         
          %%%%% removes n's when there is a mt = -1 presence
     
     for CountF = 1:length(MatE)
         NRemove2 = (MatE(:,5:8) == 1);
         MtRemove2 = (MatE(:,1:4) == 1);
     end
         
     NandMt = [NRemove2 MtRemove2];
     NandMtSum = sum(NandMt,2);
     
     Cg = 0;
     for CountG = 1:length(NandMtSum)
         if NandMtSum(CountG) < 2
             Cg = Cg + 1;
             MatF(Cg,:) = MatE(CountG,:); 
         end      
     end

         
         for CountF = 1:length(MatF)
             CountF
         k1 = MatF(CountF,1);
         k2 = MatF(CountF,2);
         l1 = MatF(CountF,3);
         l2 = MatF(CountF,4);
        
         nk1 = MatF(CountF,5);
         nk2 = MatF(CountF,6);
         nl1 = MatF(CountF,7);
         nl2 = MatF(CountF,8);
         
         mk1 = MatF(CountF,9);
         mk2 = MatF(CountF,10);
         ml1 = MatF(CountF,11);
         ml2 = MatF(CountF,12);
         
         
         
         LagTable(CountF,1) = k1;
         LagTable(CountF,2) = k2;
         LagTable(CountF,3) = l1;
         LagTable(CountF,4) = l2;
         LagTable(CountF,5) = nk1;
         LagTable(CountF,6) = nk2;
         LagTable(CountF,7) = nl1;
         LagTable(CountF,8) = nl2;
         LagTable(CountF,9) = mk1;
         LagTable(CountF,10) = mk2;
         LagTable(CountF,11) = ml1;
         LagTable(CountF,12) = ml2;
         
               
         %%%%%%%% METHOD A
%          PiTerm = sqrt(1/(factorial(nk1+abs(mk1))*factorial(nk2+abs(mk2))*factorial(nl1+abs(ml1))*factorial(nl2+abs(ml2))));
%          
%          Mt = (abs(mk1)+abs(mk2)+abs(ml1)+abs(ml2))/2;
%          fun = @(x) exp(-x).*(x.^(Mt)).* (laguerreL(nk1,abs(mk1),(x/2))) .* (laguerreL(nk2,abs(mk2),(x/2))) .* (laguerreL(nl1,abs(ml1),(x/2))) .* (laguerreL(nl2,abs(ml2),(x/2)));
%          
%          I2func = integral(fun,0,inf);
%          
%          MomTerm = 1/(2^((abs(mk1)+abs(mk2)+abs(ml1)+abs(ml2))/2));
%          
%          TotalConst = PiTerm.*I2func.*MomTerm;
%          
%          LagTable(CountF,13) = TotalConst;
         %%%%%%%%%%%%%%%%%
         %%%%%%%%% METHOD B
         MomTerm = 1/(2^((abs(mk1)+abs(mk2)+abs(ml1)+abs(ml2))/2));
    
    fun = @(x) exp(-x).*x.^((abs(mk1)+abs(mk2)+abs(ml1)+abs(ml2))/2).*...
        laguerreL(nk1,abs(mk1),(x/2)).*laguerreL(nk2,abs(mk2),(x/2)).*...
        laguerreL(nl1,abs(ml1),(x/2)).*laguerreL(nl2,abs(ml2),(x/2));
    I2func = integral(fun,0,inf);
    
    SqrtTerm = sqrt(1/(factorial(nk1 + abs(mk1))*factorial(nk2 + abs(mk2))*factorial(nl1+abs(ml1))*factorial(nl2+abs(ml2))));
    
    Const = SqrtTerm.*I2func.*MomTerm;
    
    LagTable(CountF,13) = Const;
    LagTable(CountF,14) = I2func;
    LagTable(CountF,15) = SqrtTerm;
    LagTable(CountF,16) = MomTerm;
         
         end
         
         csvwrite('LagIntLL2.csv',LagTable);
   
         
     end
              
         
     
     
     
     
     
     