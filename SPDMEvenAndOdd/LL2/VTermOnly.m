%%%%% THis method finds the allowed k values of the interaction matrix First!
tic
clc; clear all;

KPos = [ 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22; %Matlab position
    0 0 0 0 0 0 0 0 0 0  0  1  1  1  1  1  1  1  1  1  1  1 ; %n value
    -1 0 1 2 3 4 5 6 7 8  9 -1  0  1  2  3  4  5  6  7  8  9]; %Mt value

format long;
AllBras = csvread('LuisBasisLL2.csv');
%AllBras = csvread('EvenBasisLL2.csv');
AllKets = AllBras;
VMatnoATwo = zeros(length(AllBras));
aa=0;
for VKets = 1:length(AllKets)
    InitialKet = AllKets(VKets,:);
    AllowedKVal = KSelectionVMat(InitialKet); %%% FUNCTION!
    [x, y] = size(AllowedKVal);
    
    for CountA = 1:x
        
        TestKet = InitialKet;       
        k1 = AllowedKVal(CountA,1);
        k2 = AllowedKVal(CountA,2);       
        nk1 = AllowedKVal(CountA,3);
        nk2 = AllowedKVal(CountA,4);       
        mk1 = AllowedKVal(CountA,5);
        mk2 = AllowedKVal(CountA,6);       
        Opk1 = sqrt(TestKet(1,k2));
        TestKet(1,k2) = TestKet(1,k2) -1;
        Opk2 = sqrt(TestKet(1,k1)+1);
        TestKet(1,k1) = TestKet(1,k1) +1;
        TestKet;
        Op = Opk1*Opk2;
        
        for VBras = 1:length(AllBras)
            PosTester(1,1) = VMatnoATwo(9,3);
            if TestKet == AllBras(VBras,:);
                TestKet;
                CountA;
                VKets;
                VBras;
                nk1;
                nk2;
                K = (abs(mk1)+abs(mk2)+2)/2;
                Op;
                
                RootA = factorial(nk1)*factorial(nk2);
                RootB = factorial(nk1 + abs(mk1))*factorial(nk2 + abs(mk2));
                Root = sqrt(RootA/RootB);
                %fun = @(x) exp(-x).*(x.^K).* (laguerreL(nk1,abs(mk1),(x))).*(laguerreL(nk2,abs(mk2),(x)));
                if nk1 + nk2 == 0
                    I1 = factorial((abs(mk1)+abs(mk2)+2)/2);
                elseif nk1 + nk2 == 1
                    if nk1 == 1
                        MofN = abs(mk1);
                    elseif nk2 == 1
                        MofN = abs(mk2);
                    end
                    I1 = (factorial(K)*(1+MofN)) - factorial(K + 1);
                    %I1 = integral(fun,0,inf);
                   
                    
                elseif nk1 + nk2 == 2
                    I1 = factorial(K)*(1+abs(mk1)+abs(mk2)+(abs(mk1)*abs(mk2)))...
                        -((factorial(K+1))*(2+abs(mk2)+abs(mk1)))...
                        + factorial(K+2);
                    %I1 = integral(fun,0,inf);                  
                    %fun = @(x) exp(-x).*(x.^(Mt)).* (laguerreL(nk1,abs(mk1),(x/2))) .* (laguerreL(nk2,abs(mk2),(x/2))) .* (laguerreL(nl1,abs(ml1),(x/2))) .* (laguerreL(nl2,abs(ml2),(x/2)));              
                end            
                VMatnoATwo(VBras,VKets) = VMatnoATwo(VBras,VKets) + I1*Root*Op;
                if PosTester(1,1) ~= VMatnoATwo(9,3)
                    aa = aa + 1;
                    TestMat(aa,1) = VKets;
                    TestMat(aa,2) = VBras;
                end
            end
        end
        VMatnoA = 0.5.*VMatnoATwo;    
    end
end
format long
csvwrite('VMatnoA.csv',VMatnoA);
toc