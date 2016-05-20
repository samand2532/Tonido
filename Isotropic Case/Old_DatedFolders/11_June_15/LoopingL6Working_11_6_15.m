clear all;
Lmax = 6;
%%%%%%%%%%%%%%%%%% GENERATES ALL POSSIBLE M values upto Mt = 8.
StorageMatrixAll = zeros(4^Lmax,4);
StorageMatrixDelta = zeros(4^Lmax,4);
n=0;
l=0;
g=1;
f='fail';
%%%%% Generates all arrangements of m1:m4.
for m1 = 1:(Lmax+1);
    for m2 = 1:(Lmax+1);
        for m3 = 1:(Lmax+1);
            for m4 = 1:(Lmax+1);
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
%%%%%%c = cleaned matrix, excluding 0 0 0 0.
for p=1:(4^Lmax)
    mL = StorageMatrixAll(p,1)+StorageMatrixAll(p,2);
    mR = StorageMatrixAll(p,3)+StorageMatrixAll(p,4);
    if (mL==mR & mL >0 & mL <= Lmax)
        q=q+1;
        c(q,:)=StorageMatrixAll(p,:);
    end
    
end
%%%%% This section inserts 0 0 0 0 once.
AllPossibleArrangements = [0 0 0 0; c];

%%%%%%%%%%%% Big constant time%%%%%%%
%BigConst = (factorial(m1 + m2)/ sqrt(factorial(m1)*factorial(m2)*factorial(m3)*factorial(m4))) * (1/ 2^(m1+m2+1))*(g/p);
for q=1:length(AllPossibleArrangements)
    
    AllPossibleArrangements(q,5) = (factorial(AllPossibleArrangements(q,1) + AllPossibleArrangements(q,2))/ sqrt(factorial(AllPossibleArrangements(q,1))*factorial(AllPossibleArrangements(q,2))*factorial(AllPossibleArrangements(q,3))*factorial(AllPossibleArrangements(q,4)))) * (1/ 2^(AllPossibleArrangements(q,1)+AllPossibleArrangements(q,2)+1))*(g/pi);
end
%%%%%%AllPossibleArrangements Has Form |m1|m2|m3|m4|BigConstant
%                                      |  |  |  |  |

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Setting up states
InitialBras = [ 5 0 0 0  1 0 0 ;...
                4 1 0 1 0 0 0;...
                4 0 2 0   0 0 0;...
                3 2 1  0 0 0 0;...
                2 4  0 0 0 0 0];
                



InitialKets = [5 4 4 3 2 ;...
               4 1 0 1 0 ;...
               4 0 2 0 0;...
               3 2 1 0 0;...
               2 4 0 0 0;...
               0 0 0 0 0;...
               0 0 0 0 0];
    
    
    for kets=1:5
     for bras=1:5
        
        
        %InitialKets = InitialBras'
        InitialBra = InitialBras(bras,:)';
        %InitialKet = InitialBras(kets,:)';
        InitialKet = InitialKets(:,kets);
        
        %%%%%%%%%%%%%%% m ann / crea vectors
        for z=1:length(AllPossibleArrangements);
            InitialKetTest1 = 0;
            InitialKetTest2 = 0;
            Orth=0;
            Constm4Tmp=0;
            Constm3Tmp=0;
            Constm4=0;
            Constm3=0;
            Constm2=0;
            Constm1=0;
            %%%%%All variables TEST AREA!!!!!!!%%%%
            trans1=0;
            trans4=0;
            m1=0;
            m1Vect=0;
            m2=0;
            m2Vect=0;
            m3=0;
            m3Vect=0;
            m4=0;
            m4Vect=0;
            ConstComb=0;
            
            m1 = AllPossibleArrangements(z,1); m2 = AllPossibleArrangements(z,2); m3 =AllPossibleArrangements(z,3); m4 = AllPossibleArrangements(z,4);
            m1Vect = zeros(7,1); m2Vect = zeros(7,1); m3Vect = zeros(7,1); m4Vect = zeros(7,1);
            m1Vect((m1+1),1)=-1;m2Vect((m2+1),1)=-1;m3Vect((m3+1),1)=-1;m4Vect((m4+1),1)=-1;
            
            if InitialKet((m4+1),1)+m4Vect >= 0
                InitialKetTest1 = 1;
                Constm4Tmp = sqrt(InitialKet((m4+1),1));
                trans4 = InitialKet+m4Vect;
                if trans4((m3+1),1)+m3Vect >=0
                    InitialKetTest2 = 1;
                    Constm3Tmp = sqrt(trans4((m3+1),1));
                    FinalKet = trans4+m3Vect;
                else disp(f)
                end
            end
            
            if InitialKetTest1 + InitialKetTest2 == 2
                Constm4 = Constm4Tmp;
                Constm3 = Constm3Tmp;
                Constm1 = sqrt(InitialBra((m1+1),1));
                trans1 = InitialBra + m1Vect;
                Constm2 = sqrt(trans1((m2+1),1));
                FinalBra = trans1 + m2Vect;
            else
                Constm4 = 0;
                Constm3 = 0;
                Constm2 = 0;
                Constm1 = 0;
            end
            
            %%%%%%%%Orthoganality Test
            if (FinalBra(1,1) == FinalKet(1,1)) && (FinalBra(2,1) == FinalKet(2,1)) && (FinalBra(3,1) == FinalKet(3,1)) && (FinalBra(4,1) == FinalKet(4,1)) && (FinalBra(5,1) == FinalKet(5,1)) && (FinalBra(6,1) == FinalKet(6,1)) && (FinalBra(7,1) == FinalKet(7,1));
                Orth = 1;
            else Orth = 0;
            end
            %%%%%%%%If orthogonal find constants and save to matrix
            if Orth == 1;
                if (Constm4 ~=0 && Constm3 ~=0 && Constm2 ~=0 &&Constm1 ~=0)
                    ConstComb = Constm4*Constm3*Constm2*Constm1;
                    Matelement1(z,1) = ConstComb*AllPossibleArrangements(z,5);
                end
            else Matelement1(z,1)=0;
            end
            
            
            
        end
        UtermMat(bras,kets) = 0.5*sum(Matelement1);
    end
    end
    
    UtermMat
    
    
    
