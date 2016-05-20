function [  ] = PotentialLL2( LengthVect, LLVal )
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
elseif LLVal ==2
    AllBras = csvread('StatesLLL_and_LL2.csv');
end
AllKets = AllBras';

delta = [4 1; 5 1; 6 2;6 3;7 2;8 4;8 5;9 4; 10 6; 10 7;11 6;12 8;12 9;13 8;14 10;14 11;15 10;16 12;16 13;17 12;18 14;18 15;19 14;...
    1 4; 1 5;2 6;2 7; 3 6;4 8;4 9;5 8; 6 10;6 11;7 10;8 12;8 13;9 12;10 14;10 15;11 14;12 16;12 17;13 16;14 18;14 19;15 18];
%%%% Looping K to find mk1, mk+/- 2
     %%%% Looping K to find mk1, mk+/- 2
     aa=0;bb=0;
% for kk1 = 1:19
%     for kk2 = 1:19
%         mk1 = KPos(3,kk1); mk2 = KPos(3,kk2);
%         
%         if mk1 == (mk2+2)
%             aa=aa+1;
%             deltaP(aa,1) = kk1;
%             deltaP(aa,2) = kk2;
%         end
%         if mk1 ==(mk2-2)
%             bb=bb+1;
%             deltaM(bb,1) = kk1;
%             deltaM(bb,2) = kk2;
%         end
%     end
% end
% delta = [deltaP;deltaM];
   
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
                
                Mt = (abs(mk1)+abs(mk2)+2)/2;
                fun = @(x) exp(-x).*(x.^(Mt)).* (laguerreL(nk1,abs(mk1),(x))) .* (laguerreL(nk2,abs(mk2),(x)));
    
                I1func = integral(fun,0,inf);
                
                SqrtTerm = sqrt (  1/ (factorial(nk1 + abs(mk1)) * factorial(nk2 + abs(mk2))));
                
                TotalConsts = OpConst*I1func* SqrtTerm ;
                VMat(BrasV,KetsV) = VMat(BrasV,KetsV) + TotalConsts ;
            end
            
        end
    end
end
csvwrite('VMatNoA_LL2.csv',VMat);

