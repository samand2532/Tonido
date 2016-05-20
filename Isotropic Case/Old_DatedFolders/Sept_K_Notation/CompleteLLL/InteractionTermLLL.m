clc; clear all;
tic
KPos = [ 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19; %Matlab position
          0 0 1 0 1 0 1 0 1 0  1  0  1  0  1  0  1  0  1 ; %n value
         -1 0 0 1 1 2 2 3 3 4  4  5  5  6  6  7  7  8  8]; %Mt value
% k(n,m)
%           mt-1   | mt0         | mt1         | mt2         | mt3          | ...
%           (0,-1) | (0,0) (1,0) | (0,1) (1,1) | (0,2) (1,2) |  (0,3) (1,3) | ...
% Mat posit   1        2     3       4     5       6     7        8     9

g=1;
MomValue = [ -1 0 0 1 1 2 2 3 3 4 4 5 5 6 6 7 7 8 8];
NPosVect = [  0 0 1 0 1 0 1 0 1 0 1 0 1 0 1 0 1 0 1];

AllBras = csvread('StatesLLL.csv');
AllKets = AllBras';
delta = csvread('DeltaTermAllowedPositionLLL.csv');
UMatNoConst = zeros(length(AllKets));

annl1Vect = zeros(19,1);
annl2Vect = zeros(19,1);
creak1Vect = zeros(19,1);
creak2Vect = zeros(19,1);

for BrasU = 1:length(AllBras)
    for KetsU = 1:length(AllBras)
        InitialBraU = AllBras(BrasU,:);
        InitialKetU = AllKets(:,KetsU);
        
        for aa = 1:length(delta)
            annl2Const = 0;
            annl1Const = 0;
            creak2Const = 0;
            creak1Const = 0;
            
            k1 = 0; k2 = 0; l1 = 0; l2 = 0;
            annl1Vect = annl1Vect.*0;
            annl2Vect = annl2Vect.*0;
            creak1Vect = creak1Vect.*0;
            creak2Vect = creak2Vect.*0;
            
            
            k1 = delta(aa,1);
            k2 = delta(aa,2);
            l1 = delta(aa,3);
            l2 = delta(aa,4);
            annl2Vect(l2,1) = -1;
            annl1Vect(l1,1) = -1;
            creak2Vect(k2,1) = 1;
            creak1Vect(k1,1) = 1;
            
            %%%Ann / Crea Part
            %%% NEED TO STOP n>1 !!!! VERY IMPORTANT!!!!
            annl2Const = sqrt(InitialKetU(l2,1));
            l2trans = InitialKetU + annl2Vect;
            annl1Const = sqrt(l2trans(l1,1));
            l1trans = l2trans + annl1Vect;
            creak2Const = sqrt(l1trans(k2,1)+1);
            k2trans = l1trans + creak2Vect;
            creak1Const = sqrt(k2trans(k1,1)+1);
            FinalKet = k2trans + creak1Vect;
           
            
            if FinalKet == InitialBraU'
                
                OpConst = annl2Const*annl1Const*creak2Const*creak1Const;
                
                % Pi Term INC SQRT
                
                nk1 = KPos(2,k1);
                nk2 = KPos(2,k2);
                nl1 = KPos(2,l1); nl2 = KPos(2,l2);
                
                mk1 = KPos(3,k1);
                mk2 = KPos(3,k2);
                ml1 = KPos(3,l1);
                ml2 = KPos(3,l2);
                
%                 PiInt = (factorial(nk1)*factorial(nk2)*factorial(nl1)*factorial(nl2)) /...
%                     ((factorial(nk1 + abs(mk1)))*(factorial(nk2 + abs(mk2)))*(factorial(nl1 + abs(ml1)))*(factorial(nl2 + abs(ml2))));
%                 PiTerm = sqrt(PiInt);
                PiTerm = sqrt(1 / (factorial(mk1) * factorial(mk2) * factorial(ml1) * factorial(ml2)));
                %I2 = factorial(  ( abs(mk1) + abs(mk2) +abs(ml1) + abs(ml2)  ) /2    )
                %%% Momentum Term
                
                Mt = (abs(mk1)+abs(mk2)+abs(ml1)+abs(ml2))/2;
                fun = @(x) exp(-x).*(x.^(Mt)).* (laguerreL(nk1,abs(mk1),(x/2))) .* (laguerreL(nk2,abs(mk2),(x/2))) .* (laguerreL(nl1,abs(ml1),(x/2))) .* (laguerreL(nl2,abs(ml2),(x/2)));
    
                I2func = integral(fun,0,inf);
                
                Mom = 1/(2^((abs(mk1)+abs(mk2)+abs(ml1)+abs(ml2))/2));
                
                TotalConst = OpConst*PiTerm*Mom*I2func;
                UMatNoConst(BrasU,KetsU) = UMatNoConst(BrasU,KetsU) + TotalConst;
            end
            
        end
    end
end
UMatFinal = (g/(4*pi))*UMatNoConst;
csvwrite('UMatFinalLLL.csv',UMatFinal);



toc
