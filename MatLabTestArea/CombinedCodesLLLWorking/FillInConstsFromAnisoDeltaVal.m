clc; clear all;

tic
format short

g = 1;
N = 6;

KPos = [ 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19; %Matlab position
         0 0 1 0 1 0 1 0 1 0  1  0  1  0  1  0  1  0  1 ; %n value
        -1 0 0 1 1 2 2 3 3 4  4  5  5  6  6  7  7  8  8]; %Mt value

%AnisoDelta = csvread('AnisoDelta.csv');
AnisoDelta = csvread('NewAnisoDeltaNoConst.csv');

%%%% Fills the delta matrix with the constants that depend on its values
%%%% col13 = MomTerm col14 = sqrt, col15 = I2
for DeltaFill = 1:length(AnisoDelta)
    
            k1 = AnisoDelta(DeltaFill,1);
            k2 = AnisoDelta(DeltaFill,2);
                        
            nk1 = AnisoDelta(DeltaFill,3);
            nk2 = AnisoDelta(DeltaFill,4);
                   
            mk1 = AnisoDelta(DeltaFill,5);
            mk2 = AnisoDelta(DeltaFill,6);
          
    
                MtTerm = (abs(mk1)+abs(mk2)+2)/2;
                fun = @(x) exp(-x).*(x.^(MtTerm)).* (laguerreL(nk1,abs(mk1),(x))) .*...
                    (laguerreL(nk2,abs(mk2),(x)));
                
                I2func = integral(fun,0,inf);
                
                RootA = factorial(nk1)*factorial(nk2);
                RootB = factorial(nk1+abs(mk1))*factorial(nk2+abs(mk2));
                RootTerm = sqrt(RootA/RootB);
                
                
                AnisoDelta(DeltaFill,7) = RootTerm;
                AnisoDelta(DeltaFill,8) = I2func;


end


toc