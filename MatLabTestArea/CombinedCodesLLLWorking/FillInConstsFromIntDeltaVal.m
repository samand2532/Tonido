clc; clear all;

tic
format short
A = 0.00;
g = 1;
N = 6;

KPos = [ 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19; %Matlab position
         0 0 1 0 1 0 1 0 1 0  1  0  1  0  1  0  1  0  1 ; %n value
        -1 0 0 1 1 2 2 3 3 4  4  5  5  6  6  7  7  8  8]; %Mt value

IntDelta = csvread('deltaInt.csv');
%%%% Fills the delta matrix with the constants that depend on its values
%%%% col13 = MomTerm col14 = sqrt, col15 = I2
for DeltaFill = 1:length(IntDelta)
    
             k1 = IntDelta(DeltaFill,1);
            k2 = IntDelta(DeltaFill,2);
            l1 = IntDelta(DeltaFill,3);
            l2 = IntDelta(DeltaFill,4);
            
            nk1 = IntDelta(DeltaFill,5);
            nk2 = IntDelta(DeltaFill,6);
            nl1 = IntDelta(DeltaFill,7);
            nl2 = IntDelta(DeltaFill,8);
            
            mk1 = IntDelta(DeltaFill,9);
            mk2 = IntDelta(DeltaFill,10);
            ml1 = IntDelta(DeltaFill,11);
            ml2 = IntDelta(DeltaFill,12);
    
    
    
                halfMt = (abs(mk1)+abs(mk2)+abs(ml1)+abs(ml2))/2;
                Momterm = 1/(2^halfMt);
                
                fun = @(x) exp(-x).*(x.^(halfMt)).* (laguerreL(nk1,abs(mk1),(x/2))) .*...
                    (laguerreL(nk2,abs(mk2),(x/2))) .* (laguerreL(nl1,abs(ml1),(x/2))) .*...
                    (laguerreL(nl2,abs(ml2),(x/2)));
                
                I2func = integral(fun,0,inf);
                
                RootA = factorial(nk1)*factorial(nk2)*factorial(nl1)*factorial(nl2);
                RootB = factorial(nk1+abs(mk1))*factorial(nk2+abs(mk2))*...
                    factorial(nl1+abs(ml1))*factorial(nl2+abs(ml2));
                RootTerm = sqrt(RootA/RootB);
                
                IntDelta(DeltaFill,13) = Momterm;
                IntDelta(DeltaFill,14) = RootTerm;
                IntDelta(DeltaFill,15) = round(I2func);


end


toc