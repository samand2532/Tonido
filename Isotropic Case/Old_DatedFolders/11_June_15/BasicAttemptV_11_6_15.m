%%%% VtermTry3 - 11/6/15
clear all; clc
A=1;

AllBras = [ 6 0 0 0 0 0 0 9 9 9 9;...
            4 2 0 0 0 0 0 9 9 9 9;...
            5 0 1 0 0 0 0 9 9 9 9;...
            5 0 0 0 1 0 0 9 9 9 9;...
            4 1 0 1 0 0 0 9 9 9 9;...
            4 0 2 0 0 0 0 9 9 9 9;...
            3 2 1 0 0 0 0 9 9 9 9;...
            2 4 0 0 0 0 0 9 9 9 9];
        
AllKets = [ 6 4 5 5 4 4 3 2;...
            0 2 0 0 1 0 2 4;...
            0 0 1 0 0 2 1 0;...        
            0 0 0 0 1 0 0 0;...     
            0 0 0 1 0 0 0 0;...    
            0 0 0 0 0 0 0 0;...      
            0 0 0 0 0 0 0 0;...
            9 9 9 9 9 9 9 9;...
            9 9 9 9 9 9 9 9;...
            9 9 9 9 9 9 9 9;...
            9 9 9 9 9 9 9 9;...
            9 9 9 9 9 9 9 9];
            
VFinal = zeros(8);


for bra = 3
    for ket = 3
        
        InitialKet = AllKets(:,ket);
        InitialBra = AllBras(bra,:)';
        %%%%%%%%m0
        P2m0ConstInitial = sqrt(2) * sqrt((InitialKet(1,1)+1)) * sqrt(InitialKet(3,1));
        P2m0Ket = InitialKet + [1 0 0 0 0 0 0]' + [0 0 -1 0 0 0 0]';
        
        %%%%%m1
        P2m1ConstInitial = sqrt(6) * sqrt((InitialKet(2,1)+1)) * sqrt(InitialKet(4,1));
        P2m1Ket = InitialKet + [0 1 0 0 0 0 0]' + [0 0 0 -1 0 0 0]';
        
        %%%%m2
        P1m2ConstInitial = sqrt(2) * sqrt((InitialKet(3,1)+1)) * sqrt(InitialKet(1,1));
        P1m2Ket = InitialKet + [0 0 1 0 0 0 0]' + [-1 0 0 0 0 0 0]';
        P2m2ConstInitial = sqrt(12) * sqrt((InitialKet(3,1)+1)) * sqrt(InitialKet(5,1));
        P2m2Ket = InitialKet + [0 0 1 0 0 0 0]' + [0 0 0 -1 0 0 0]';
        
        %%%%%m3
        P1m3ConstInitial = sqrt(6) * sqrt((InitialKet(4,1)+1)) * sqrt(InitialKet(2,1));
        P1m3Ket = InitialKet + [0 0 0 1 0 0 0]' + [0 -1 0 0 0 0 0]';
        P2m3ConstInitial = sqrt(20) * sqrt((InitialKet(4,1)+1)) * sqrt(InitialKet(6,1));
        P2m3Ket = InitialKet + [0 0 0 1 0 0 0]' + [0 0 0 0 0 -1 0]';
        
        %%%%%m4
        P1m4ConstInitial = sqrt(12) * sqrt((InitialKet(5,1)+1)) * sqrt(InitialKet(3,1));
        P1m4Ket = InitialKet + [0 0 0 0 1 0 0]' + [0 0 0 0 -1 0 0]';
        P2m4ConstInitial = sqrt(30) * sqrt((InitialKet(5,1)+1)) * sqrt(InitialKet(7,1));
        P2m4Ket = InitialKet + [0 0 0 0 1 0 0]' + [0 0 0 0 0 0 -1]';
        
        %%%%%%m5
        P1m5ConstInitial = sqrt(20) * sqrt((InitialKet(6,1)+1)) * sqrt(InitialKet(4,1));
        P1m5Ket = InitialKet + [0 0 0 0 0 1 0]' + [0 0 0 -1 0 0 0]';
        %P2m5ConstInitial = sqrt(42) * sqrt((InitialKet(6,1)+1)) * sqrt(InitialKet(8,1));
        %P2m5Ket = InitialKet + [0 0 0 0 0 1 0]' + [0 0 0 0 0 0 0]';
        
        
        
        
        [P2m0Ortho] = Orthogonal_11_6_15 (P2m0Ket, InitialBra);
        [P2m1Ortho] = Orthogonal_11_6_15 (P2m0Ket, InitialBra);
        [P1m2Ortho] = Orthogonal_11_6_15 (P2m0Ket, InitialBra);
        [P2m2Ortho] = Orthogonal_11_6_15 (P2m0Ket, InitialBra);
        [P1m3Ortho] = Orthogonal_11_6_15 (P2m0Ket, InitialBra);
        [P2m3Ortho] = Orthogonal_11_6_15 (P2m0Ket, InitialBra);
        [P1m4Ortho] = Orthogonal_11_6_15 (P2m0Ket, InitialBra);
        [P2m4Ortho] = Orthogonal_11_6_15 (P2m0Ket, InitialBra);
        [P1m5Ortho] = Orthogonal_11_6_15 (P2m0Ket, InitialBra);
        [P2m5Ortho] = Orthogonal_11_6_15 (P2m0Ket, InitialBra);
        if P2m0Ortho ==1
            VFinal(bra,ket) = VFinal(bra,ket) + P2m0ConstInitial;
        end   
        if P2m1Ortho ==1
            VFinal(bra,ket) = VFinal(bra,ket) + P2m1ConstInitial;
        end 
        if P1m2Ortho ==1
            VFinal(bra,ket) = VFinal(bra,ket) + P1m2ConstInitial;
        end  
        if P2m2Ortho ==1
            VFinal(bra,ket) = VFinal(bra,ket) + P2m2ConstInitial;
        end  
        if P1m3Ortho ==1
            VFinal(bra,ket) = VFinal(bra,ket) + P1m3ConstInitial;
        end  
        if P2m3Ortho ==1
            VFinal(bra,ket) = VFinal(bra,ket) + P2m3ConstInitial;
        end  
        if P1m4Ortho ==1
            VFinal(bra,ket) = VFinal(bra,ket) + P1m4ConstInitial;
        end  
        if P2m4Ortho ==1
            VFinal(bra,ket) = VFinal(bra,ket) + P2m4ConstInitial;
        end  
        if P1m5Ortho ==1
            VFinal(bra,ket) = VFinal(bra,ket) + P1m5ConstInitial;
        end  
        %if P2m5Ortho ==1
        %    VFinal(bra,ket) = VFinal(bra,ket) + P2m5ConstInitial;
        %end  
        
    end
end
VFinal

