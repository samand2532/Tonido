clear all;
InitialBras = [5 0 0 0 0 0 1 ;...
    4 1 0 0 0 1 0 ;...
    4 0 1 0 1 0 0 ;...
    3 2 0 0 1 0 0 ;...
    4 0 0 2 0 0 0 ;...
    3 1 1 1 0 0 0 ;...
    2 3 0 1 0 0 0 ;...
    3 0 3 0 0 0 0 ;...
    2 2 2 0 0 0 0 ;...
    1 4 1 0 0 0 0 ;...
    0 6 0 0 0 0 0 ];



InitialKets = [5 4 4 3 4 3 2 3 2 1 0;...
    0 1 0 2 0 1 3 0 2 4 6;...
    0 0 1 0 0 1 0 3 2 1 0;...
    0 0 0 0 2 1 1 0 0 0 0;...
    0 0 1 1 0 0 0 0 0 0 0;...
    0 1 0 0 0 0 0 0 0 0 0;...
    1 0 0 0 0 0 0 0 0 0 0];
mloopmat = zeros(7,1);
CreaVect = zeros(7,1);
AnnaVect = zeros(7,1);
P1CreaVect = CreaVect; P2CreaVect = CreaVect;
P1AnnaVect = AnnaVect; P2AnnaVect = AnnaVect;
%%% For L = 6  mat = 11x11
Vmat = zeros(11);
mfinal=zeros(11);

% m cannot go < 1 or >7
%REMEMBER THAT HAVE TO START AT 1!!!
% p1 = root(m(m-1))*adagm*acrea(m-2)
% p2 = root((m+1)(m+2))*adagm * acrea (m+2)


for ket = 1:11
    for bra =1:1
        
        
        InitialBra = InitialBras(bra,:)';
        InitialKet = InitialKets(:,ket);
        P1=0;
        P2=0;
        
        
        
        for m = 4:4 % from 1 to 7 - FOR MATLAB - ann / crea changed
            if (m-1)>1
                P1CreaVect(m-1,1) = 1;
            else
                P1CreaVect = InitialKet.*9; %%%% .* 9 to ensure non-orthogonal
            end
            
            if (m-3)>0
                P1AnnaVect(m-3,1) = -1;
            else
                P1AnnaVect = InitialKet.*9; %%%% .* 9 to ensure non-orthogonal
            end
            if (m-1)>1
                P2CreaVect(m-1,1) = 1;
            else
                P2CreaVect = InitialKet.*9; %%%% .* 9 to ensure non-orthogonal
            end
            P2AnnaVect(m+1,1) = -1;
            
            %%%%%P1
            if (m-3) >0
                P1=sqrt((m-1)*(m-2))*sqrt(InitialKet((m-1),1))*sqrt(InitialKet((m-3),1));
                P1trans = InitialKet + P1CreaVect + P1AnnaVect;
                
            else
                P1=0;
                P1trans = InitialKet .*9;
            end
            [P1Orthog] = Orthogonal_6_6_15 (P1trans',InitialBra')
            if P1Orthog ==1;
                P1Final = P1;
            else P1Final = 0
            end
            
            %%%%%%P2
            if ((m-1)>0) & ((m+1)<8)
                P2 = sqrt((m*(m-1)))*sqrt(InitialKet((m-1),1))*sqrt(InitialKet((m+1),1));
                P2trans = InitialKet + P1CreaVect + P1AnnaVect;
            else
                P2=0;
                P2trans = InitialKet.*9;
            end
            [P2Orthog] = Orthogonal_6_6_15 (P2trans',InitialBra')
            if P2Orthog ==1;
                P2Final = P2;
            else P2Final = 0
            end
            
            mloopmat(m,1) = P1Final+P2Final;
            mfinal(bra,ket) = sum(mloopmat);
        end
    end
end






P1
P2
P1trans
P2trans

mfinal



