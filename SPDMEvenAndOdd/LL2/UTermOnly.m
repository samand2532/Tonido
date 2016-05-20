%%%%% THis method finds the allowed k values of the interaction matrix First!
tic
clc; clear all;

KPos = [ 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22; %Matlab position
    0 0 0 0 0 0 0 0 0 0  0  1  1  1  1  1  1  1  1  1  1  1 ; %n value
    -1 0 1 2 3 4 5 6 7 8  9 -1  0  1  2  3  4  5  6  7  8  9]; %Mt value


AllBras = csvread('LuisBasisLL2.csv');
%AllBras = csvread('EvenBasisLL2.csv');
AllKets = AllBras;
UMatNoConst = zeros(length(AllBras));

for UKets = 1:length(AllKets)
    InitialKet = AllKets(UKets,:);
    deltaFinal = KSelectionUMat(InitialKet); %%% FUNCTION!
    [x, y] = size(deltaFinal);
    
    for CountA = 1:x
        
        
        TestKet = InitialKet;
        InitialKet;
        
        k1 = deltaFinal(CountA,1);
        k2 = deltaFinal(CountA,2);
        l1 = deltaFinal(CountA,3);
        l2 = deltaFinal(CountA,4);
        
        nk1 = deltaFinal(CountA,5);
        nk2 = deltaFinal(CountA,6);
        nl1 = deltaFinal(CountA,7);
        nl2 = deltaFinal(CountA,8);
        
        mk1 = deltaFinal(CountA,9);
        mk2 = deltaFinal(CountA,10);
        ml1 = deltaFinal(CountA,11);
        ml2 = deltaFinal(CountA,12);
        
        OpK1 = sqrt(TestKet(1,k1));
        TestKet(1,k1) = TestKet(1,k1)-1;
        OpK2 = sqrt(TestKet(1,k2));
        TestKet(1,k2) = TestKet(1,k2)-1;
        TestKet(1,l1) = TestKet(1,l1)+1;
        OpL1 = sqrt(TestKet(1,l1));
        TestKet(1,l2) = TestKet(1,l2)+1;
        OpL2 = sqrt(TestKet(1,l2));
        TestKet;
        Op=OpK1*OpK2*OpL1*OpL2;
        
        for UBras = 1:length(AllBras)
            if TestKet == AllBras(UBras,:)
                CountA;
                deltaFinal(CountA,1:4);
                
                
                absMt = abs(mk1)+abs(mk2)+abs(ml1)+abs(ml2);
                halfMt = absMt/2;
                MomTerm= 1/(2^halfMt);
                %%%Root Term
                RootA = factorial(nk1)*factorial(nk2)*factorial(nl1)*factorial(nl2);
                RootB = (factorial(nk1+abs(mk1)))*(factorial(nk2+abs(mk2)))*(factorial(nl1+abs(ml1)))*(factorial(nl2+abs(ml2)));
                RootTerm = sqrt(RootA/RootB);
               
                NRange = deltaFinal(CountA,5:8);
               
                [nPos] = find(NRange == 1);
                MRange = deltaFinal(CountA,9:12);
                
                if nk1+nk2+nl1+nl2 == 0
                   % disp('n=0');
                    I2 = factorial(halfMt);
                elseif nk1+nk2+nl1+nl2 == 1
                    %disp('n=1');
                    n1Pos = nPos(1,1);
                    I2 = (factorial(halfMt))*(1+abs(MRange(1,n1Pos))) - 0.5*factorial(halfMt+1);
                elseif nk1+nk2+nl1+nl2 == 2
                   % disp('n=2');
                    n1Pos = nPos(1,1);
                    n2Pos = nPos(1,2);
                    I2 = (factorial(0.5*absMt))*...
                        ( (1+abs(MRange(1,n1Pos)))*(1+abs(MRange(1,n2Pos))) - 0.5*(2+abs(MRange(1,n1Pos))+abs(MRange(1,n2Pos)))*(1+0.5*absMt)+...
                        0.25*(1+0.5*absMt)*(2+0.5*absMt));
                    
                    
                    %I2 = (factorial(0.5*absMt))*(1 + abs(MRange(1,n1Pos))+abs(MRange(1,n2Pos))+...
                    %    (abs(MRange(1,n1Pos))*abs(MRange(1,n2Pos)))) - ...
                    %    factorial(0.5*absMt+1)*(1 + 0.5*abs(MRange(1,n1Pos))*0.25*abs(MRange(1,n2Pos)))+...
                    %    0.25*factorial((0.5*absMt)+2);
                end

       
                UMatNoConst(UKets,UBras) = UMatNoConst(UKets,UBras) + (MomTerm * Op * RootTerm * I2);
            end
        end
    end
end

UMat = 1/(4*pi) .* UMatNoConst;
csvwrite('MyUMat.csv',UMat);





toc