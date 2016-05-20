function [  ] = PotentialLLL( LengthVect, LLVal )
%UNTITLED6 Summary of this function goes here
%   Detailed explanation goes here
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


if LLVal == 1
AllBras = csvread('StatesLLL.csv');
elseif LLVavl ==2
    AllBras = csvread('StatesLLL_and_LL2.csv');
end
AllKets = AllBras';

delta = [2 6; 4 8; 4 6; 6 10; 8 12; 10 14; 12 16; 14 18; 6 2; 8 4; 10 6; 12 8; 14 10; 16 12; 18 14];
VMat = zeros(length(AllKets));


creak1Vect = zeros(LengthVect,1);
annk2Vect  = zeros(LengthVect,1);
for BrasV = 1:length(AllBras)
    for KetsV = 1:length(AllBras)
        InitialBraV = AllBras(BrasV,:);
        InitialKetV = AllKets(:,KetsV);
        
        for aa = 1:length(delta)
            k1 = 0; k2 = 0;
            creak1Vect = creak1Vect .* 0;
            annk2Vect = annk2Vect .* 0;
            k1 = delta(aa,1);
            k2 = delta(aa,2);
                        
            annk2Vect(k2,1) = -1;
            creak1Vect(k1,1) = 1;
            
            annk2Const = sqrt(InitialKetV(k2,1));
            k2trans = InitialKetV + annk2Vect;
            creak1Const = sqrt(k2trans(k1,1)+1);
            FinalKet = k2trans + creak1Vect;
            
            if FinalKet == InitialBraV';
                
                OpConst = annk2Const * creak1Const;
                
                nk1 = KPos(2,k1); nk2 = KPos(2,k2);
                mk1 = KPos(3,k1); mk2 = KPos(3,k2);
                
                I1 = factorial( (abs(mk1) + abs(mk2) + 2)/2 );
                SqrtTerm = sqrt (  1/ (factorial(nk1 + abs(mk1)) * factorial(nk2 + abs(mk2))));
                
                TotalConsts = OpConst*I1* SqrtTerm ;
                VMat(BrasV,KetsV) = VMat(BrasV,KetsV) + TotalConsts ;
            end
            
        end
    end
end
csvwrite('VMatNoA_LLL.csv',VMat);

