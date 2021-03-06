%%%% Vterm
clear all; clc
A=1;
         
InitialBras = [ 6 0 0 0 0;...
                4 2 0 0 0 ;...
                5 0 1 0 0 ;...
                5 0 0 0 1 ;...
                4 1 0 1 0 ;...
                4 0 2 0 0 ;...
                3 2 1 0 0 ;...
                2 4 0 0 0 ];

            
InitialKets = [ 6 4 5 5 4 4 3 2;...
                0 2 0 0 1 0 2 4;...
                0 0 1 0 0 2 1 0;...
                0 0 0 0 1 0 0 0;...
                0 0 0 1 0 0 0 0];

      

P1FinalKet = zeros((length(InitialBras(1,:))),1);

Pmmatrix=zeros(12,3);
%%%%%Pmmatrix Formated as - |P1|P2|(P1+P2) P values ONLY IF ORTHOGONAL
for bra = 1:1
    for ket = 1:1
        
        InitialKet = InitialKets(:,ket);
        InitialBra = InitialBras(bra,:)';
        
        for n=1:5
            P2ann = zeros(length(InitialKet),1);
            P1ann = zeros(length(InitialKet),1);
            crea = zeros(length(InitialKet),1);
            m=n-1;
            crea(n,1) = 1;
            P2ann(n+2,1) = -1;
            
            
            %%%%%%%%% P2
            if n+2 < 5
            P2 = sqrt((m+1)*(m+2))*sqrt(InitialKet(n,1))*sqrt(InitialKet((n+2),1));
            P2FinalKet = InitialKet + P2ann + crea;
            [P2Orth] = Orthogonal_7_6_15 (InitialBra, P2FinalKet)
            if P2Orth == 1
                Pmmatrix(n,2) = P2
            else
                Pmmatrix(n,2) = 0;
                
            end
            end
            %%%%% P1
            if (n-2) > 0
                P1ann(n-2,1) = -1;
                P1 = sqrt(m*(m-1))*sqrt(InitialKet(n,1))*sqrt(InitialKet((n-2),1));
                P1FinalKet = InitialKet + P1ann + crea
                [P1Orth] = Orthogonal_7_6_15 (InitialBra, P1FinalKet)
                if P1Orth == 1
                Pmmatrix(n,1) = P1
                else
                Pmmatrix(n,1) = 0; 
            
                
                end
                
            end
            
            
            
        end
        sum(sum(Pmmatrix))*A
        
        
        
        
    end
end
%P2
%P1
%InitialKet
%P1FinalKet


% %%%%%%%%%%%%    P1
%             if (m-2) > 0
%                 P1 = (sqrt(m*(m-1)))*sqrt(InitialKet(m,1))*sqrt(InitialKet((m-2),1));
%                 P1FinalKet = InitialKet + crea + ann;
%                 P1mmatrix(n,1) = P1
%     
%             end   