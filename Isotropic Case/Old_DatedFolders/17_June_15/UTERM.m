%NEEDS WORK SECTION!!!!!!!!!!!!!!!!!!!!!!!!!
%
%
%
%
%

%%%%%%%N TERM WORKING!!!! L TERM WORKING !!!!!!!!!!!!
clear all; clc
%%%%%%%%%%%%%%%Changing Constants
Lmax=8;
Omega = 1;
g=1;
N=6;
Leyes=zeros(39);
MomValue = [ 0 1 2 3 4 5 6 7 8 9 10]';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
AllBras = [  6     0     0     0     0     0     0     0     0     0     0;...
    5     0     1     0     0     0     0     0     0     0     0;...
    4     2     0     0     0     0     0     0     0     0     0;...
    5     0     0     0     1     0     0     0     0     0     0;...
    4     1     0     1     0     0     0     0     0     0     0;...
    4     0     2     0     0     0     0     0     0     0     0;...
    3     2     1     0     0     0     0     0     0     0     0;...
    2     4     0     0     0     0     0     0     0     0     0;...
    5     0     0     0     0     0     1     0     0     0     0;...
    4     1     0     0     0     1     0     0     0     0     0;...
    4     0     1     0     1     0     0     0     0     0     0;...
    3     2     0     0     1     0     0     0     0     0     0;...
    4     0     0     2     0     0     0     0     0     0     0;...
    3     1     1     1     0     0     0     0     0     0     0;...
    2     3     0     1     0     0     0     0     0     0     0;...
    3     0     3     0     0     0     0     0     0     0     0;...
    2     2     2     0     0     0     0     0     0     0     0;...
    1     4     1     0     0     0     0     0     0     0     0;...
    0     6     0     0     0     0     0     0     0     0     0;...
    5     0     0     0     0     0     0     0     1     0     0;...
    4     1     0     0     0     0     0     1     0     0     0;...
    4     0     1     0     0     0     1     0     0     0     0;...
    3     2     0     0     0     0     1     0     0     0     0;...
    4     0     0     1     0     1     0     0     0     0     0;...
    3     1     1     0     0     1     0     0     0     0     0;...
    2     3     0     0     0     1     0     0     0     0     0;...
    4     0     0     0     2     0     0     0     0     0     0;...
    3     1     0     1     1     0     0     0     0     0     0;...
    3     0     2     0     1     0     0     0     0     0     0;...
    2     2     1     0     1     0     0     0     0     0     0;...
    1     4     0     0     1     0     0     0     0     0     0;...
    3     0     1     2     0     0     0     0     0     0     0;...
    2     2     0     2     0     0     0     0     0     0     0;...
    2     1     2     1     0     0     0     0     0     0     0;...
    1     3     1     1     0     0     0     0     0     0     0;...
    0     5     0     1     0     0     0     0     0     0     0;...
    2     0     4     0     0     0     0     0     0     0     0;...
    1     2     3     0     0     0     0     0     0     0     0;...
    0     4     2     0     0     0     0     0     0     0     0];

AllKets =  AllBras';

NmatFinal = zeros(39);
NmatTemp = zeros(11,1);
annVect  = zeros(11,1); creaVect = annVect;

%%%%%%%%%%%%%START U TERM%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%% GENERATES ALL POSSIBLE m arrangements upto Mt = 8.
StorageMatrixAll = zeros(4^Lmax,4);
StorageMatrixDelta = zeros(4^Lmax,4);
n=0;
l=0;
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

for KetsU=1:39
    for BrasU=1:39

        InitialBraU = AllBras(BrasU,:)';
        InitialKetU = AllKets(:,KetsU);
        
        %%%%%%%%%%%%%%% m ann / crea vectors
        for z=1:length(AllPossibleArrangements);
            Constm4Tmp=0;
            Constm3Tmp=0;
            Constm4=0;
            Constm3=0;
            Constm2=0;
            Constm1=0;
            
            m1 = AllPossibleArrangements(z,1); m2 = AllPossibleArrangements(z,2); m3 =AllPossibleArrangements(z,3); m4 = AllPossibleArrangements(z,4);
            m1Vect = zeros(11,1); m2Vect = zeros(11,1); m3Vect = zeros(11,1); m4Vect = zeros(11,1);
            m1Vect((m1+1),1)=1;m2Vect((m2+1),1)=1;m3Vect((m3+1),1)=-1;m4Vect((m4+1),1)=-1;
            
            m4Const = sqrt(InitialKetU((m4+1),1));
            m4trans = InitialKetU + m4Vect;
            m3Const = sqrt(m4trans((m3+1),1));
            m3trans = m4trans + m3Vect;
            m2Const = sqrt(m3trans((m2+1),1)+1);
            m2trans = m3trans + m2Vect;
            m1Const = sqrt(m2trans((m1+1),1)+1);
            FinalKetU = m2trans + m1Vect;
            
            %%%%%%%%Orthoganality Test
            [OrthU] = Orthogonal(InitialBraU,FinalKetU);
            %%%%%%%%If orthogonal find constants and save to matrix
            
            ConstComb = OrthU*m4Const*m3Const*m2Const*m1Const;
            Matelement1(z,1) = ConstComb*AllPossibleArrangements(z,5);
            
        end
        UMatFinal(BrasU,KetsU) = 0.5*sum(Matelement1);
    end
end
disp('U done')
UMatFinal

