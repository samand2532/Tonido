clear all;

%%%%%%%%%%%%%%%%%% GENERATES ALL POSSIBLE M values upto Mt = 8.
StorageMatrixAll = zeros(1024,4);
StorageMatrixDelta = zeros(1024,4);
n=0;
l=0;
g=1;
f='fail';
%%%%% Generates all arrangements of m1:m4.
for m1 = 1:5;
    for m2 = 1:5;
        for m3 = 1:5;
            for m4 = 1:5;
                n = n+1;
                
                StorageMatrixAll(n,1)=m1-1;
                StorageMatrixAll(n,2)=m2-1;
                StorageMatrixAll(n,3)=m3-1;
                StorageMatrixAll(n,4)=m4-1;
                
            end
        end
    end
end
q=0;
%%%%% This section conserves momentum and removes 0000 multiple times.
for p=1:1024
    mL = StorageMatrixAll(p,1)+StorageMatrixAll(p,2);
    mR = StorageMatrixAll(p,3)+StorageMatrixAll(p,4);
    if (mL==mR & mL >0 & mL <= 4)
        q=q+1;
        c(q,:)=StorageMatrixAll(p,:);
    end
    
end
%%%%% This section inserts 0 0 0 0 once.
AllPossibleArrangements = [0 0 0 0; c];

%%%%%%%%%%%% Big constant time%%%%%%%
%BigConst = (factorial(m1 + m2)/ sqrt(factorial(m1)*factorial(m2)*factorial(m3)*factorial(m4))) * (1/ 2^(m1+m2+1))*(g/p);
for q=1:55
    
    AllPossibleArrangements(q,5) = (factorial(AllPossibleArrangements(q,1) + AllPossibleArrangements(q,2))/ sqrt(factorial(AllPossibleArrangements(q,1))*factorial(AllPossibleArrangements(q,2))*factorial(AllPossibleArrangements(q,3))*factorial(AllPossibleArrangements(q,4)))) * (1/ 2^(AllPossibleArrangements(q,1)+AllPossibleArrangements(q,2)+1))*(g/pi);
end
%%%%%%AllPossibleArrangements Has Form |m1|m2|m3|m4|BigConstant
%              |  |  |  |  |

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Setting up states
%InitialBras = [5 0 0 0 1 0;...
%               4 1 0 1 0 0;...
%              4 0 2 0 0 0;...
%             3 2 1 0 0 0;...
%            2 4 0 0 0 0;]
%InitialKets = InitialBras'
InitialBra = [ 6 0 0 0 0 ]';
InitialKet = [ 6 0 0 0 0 ]';

%%%%%%%%%%%%%%% m ann / crea vectors
for z=1:length(AllPossibleArrangements);
    Orth=0;
    m1 = AllPossibleArrangements(z,1); m2 = AllPossibleArrangements(z,2); m3 =AllPossibleArrangements(z,3); m4 = AllPossibleArrangements(z,4);
    m1Vect = zeros(5,1); m2Vect = zeros(5,1); m3Vect = zeros(5,1); m4Vect = zeros(5,1);
    m1Vect((m1+1),1)=-1;m2Vect((m2+1),1)=-1;m3Vect((m3+1),1)=-1;m4Vect((m4+1),1)=-1;
    
    if InitialKet((m4+1),1)+m4Vect >= 0 
        
        Constm4 = sqrt(InitialKet((m4+1),1));
        trans4 = InitialKet+m4Vect;
        if trans4((m3+1),1)+m3Vect >=0
        Constm3 = sqrt(trans4((m3+1),1));
        FinalKet = trans4+m3Vect;
        else disp(f)
        end
        Constm4
    end
    
    
    Constm1 = sqrt(InitialBra((m1+1),1));
    trans1 = InitialBra + m1Vect;
    Constm2 = sqrt(trans1((m2+1),1));
    FinalBra = trans1 + m2Vect;
    
    if (FinalBra(1,1) == FinalKet(1,1)) & (FinalBra(2,1) == FinalKet(2,1)) & (FinalBra(3,1) == FinalKet(3,1)) & (FinalBra(4,1) == FinalKet(4,1));
        Orth = 1;
    else Orth = 0;
    end
    
    if Orth == 1;
        if (Constm4 ~=0 & Constm3 ~=0 & Constm2 ~=0 &Constm1 ~=0) 
        ConstComb = Constm4*Constm3*Constm2*Constm1;
        Matelement1(z,1) = ConstComb*AllPossibleArrangements(z,5);
        end
    else
    end
    
    
    
end
Uterm = 0.5*sum(Matelement1)





