function [ I2func ] = I2Func( k1,k2,l1,l2 )
%UNTITLED9 Summary of this function goes here
%   Detailed explanation goes here

KPos = [ 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19; %Matlab position
          0 0 1 0 1 0 1 0 1 0  1  0  1  0  1  0  1  0  1 ; %n value
         -1 0 0 1 1 2 2 3 3 4  4  5  5  6  6  7  7  8  8]; %Mt value

  
    
    nk1 = KPos(2,k1); nk2 = KPos(2,k2);
    nl1 = KPos(2,l1); nl2 = KPos(2,l2);
    
    mk1 = KPos(3,k1);
    mk2 = KPos(3,k2);
    ml1 = KPos(3,l1);
    ml2 = KPos(3,l2);
    Mt = (abs(mk1)+abs(mk2)+abs(ml1)+abs(ml2))/2;
    fun = @(x) exp(-x).*(x.^(Mt)).* (laguerreL(nk1,abs(mk1)).*(x/2)) .* (laguerreL(nk2,abs(mk2)).*(x/2)) .* (laguerreL(nl1,abs(ml1)).*(x/2)) .* (laguerreL(nl2,abs(ml2)).*(x/2));
    
    I2func = integral(fun,0,inf)
    
    
    
    



end

