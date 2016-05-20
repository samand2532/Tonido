%% Vterm
AllBras = csvread('AllBras.csv');
AllKets =  AllBras';
VFinal = zeros(6);


for bra = 1:6
    for ket = 1:6
        
        InitialKet = AllKets(:,ket);
        InitialBra = AllBras(bra,:)';
        
        % m0
        P2m0Const = sqrt(2) * sqrt(InitialKet(1,1)+1) * sqrt(InitialKet(3,1));
        P2m0Ket = InitialKet + [1 0 0 0 0 0 0 0 0 0 0]' + [0 0 -1 0 0 0 0 0 0 0 0]';
        
        % m1
        P2m1Const = sqrt(6) * sqrt(InitialKet(2,1)+1) * sqrt(InitialKet(4,1));
        P2m1Ket = InitialKet + [0 1 0 0 0 0 0 0 0 0 0]' + [0 0 0 -1 0 0 0 0 0 0 0]';
        
        % m2
        P1m2Const = sqrt(2) * sqrt(InitialKet(3,1)+1) * sqrt(InitialKet(1,1));
        P1m2Ket = InitialKet + [0 0 1 0 0 0 0 0 0 0 0]' + [-1 0 0 0 0 0 0 0 0 0 0]';
        P2m2Const = sqrt(12) * sqrt(InitialKet(3,1)+1) * sqrt(InitialKet(5,1));
        P2m2Ket = InitialKet + [0 0 1 0 0 0 0 0 0 0 0]' + [0 0 0 0 -1 0 0 0 0 0 0]';
        %m3
        P1m3Const = sqrt(6) * sqrt(InitialKet(4,1)+1) * sqrt(InitialKet(2,1));
        P1m3Ket = InitialKet + [0 0 0 1 0 0 0 0 0 0 0]' + [0 -1 0 0 0 0 0 0 0 0 0]';
        P2m3Const = sqrt(20) * sqrt(InitialKet(4,1)+1) * sqrt(InitialKet(6,1));
        P2m3Ket = InitialKet + [0 0 0 1 0 0 0 0 0 0 0]' + [0 0 0 0 0 -1 0 0 0 0 0]';
        %m4
        P1m4Const = sqrt(12) * sqrt(InitialKet(5,1)+1) * sqrt(InitialKet(3,1));
        P1m4Ket = InitialKet + [0 0 0 0 1 0 0 0 0 0 0]' + [0 0 -1 0 0 0 0 0 0 0 0]';
        P2m4Const = sqrt(30) * sqrt(InitialKet(5,1)+1) * sqrt(InitialKet(7,1));
        P2m4Ket = InitialKet + [0 0 0 0 1 0 0 0 0 0 0]' + [0 0 0 0 0 0 -1 0 0 0 0]';
        %m5
        P1m5Const = sqrt(20) * sqrt(InitialKet(6,1)+1) * sqrt(InitialKet(4,1));
        P1m5Ket = InitialKet + [0 0 0 0 0 1 0 0 0 0 0]' + [0 0 0 -1 0 0 0 0 0 0 0]';
        P2m5Const = sqrt(42) * sqrt(InitialKet(6,1)+1) * sqrt(InitialKet(8,1));
        P2m5Ket = InitialKet + [0 0 0 0 0 1 0 0 0 0 0]' + [0 0 0 0 0 0 0 -1 0 0 0]';
        %m6
        P1m6Const = sqrt(30) * sqrt(InitialKet(7,1)+1) * sqrt(InitialKet(5,1));
        P1m6Ket = InitialKet + [0 0 0 0 0 0 1 0 0 0 0]' + [0 0 0 0 -1 0 0 0 0 0 0]';
        P2m6Const = sqrt(56) * sqrt(InitialKet(7,1)+1) * sqrt(InitialKet(9,1));
        P2m6Ket = InitialKet + [0 0 0 0 0 0 1 0 0 0 0]' + [0 0 0 0 0 0 0 0 -1 0 0]';
        %m7
        P1m7Const = sqrt(42) * sqrt(InitialKet(8,1)+1) * sqrt(InitialKet(6,1));
        P1m7Ket = InitialKet + [0 0 0 0 0 0 0 1 0 0 0]' + [0 0 0 0 0 -1 0 0 0 0 0]';
        P2m7Const = sqrt(72) * sqrt(InitialKet(8,1)+1) * sqrt(InitialKet(10,1));
        P2m7Ket = InitialKet + [0 0 0 0 0 0 0 1 0 0 0]' + [0 0 0 0 0 0 0 0 0 -1 0]';
        %m8
        P1m8Const = sqrt(56) * sqrt(InitialKet(9,1)+1) * sqrt(InitialKet(7,1));
        P1m8Ket = InitialKet + [0 0 0 0 0 0 0 0 1 0 0]' + [0 0 0 0 0 0 -1 0 0 0 0]';
        P2m8Const = sqrt(90) * sqrt(InitialKet(9,1)+1) * sqrt(InitialKet(11,1));
        P2m8Ket = InitialKet + [0 0 0 0 0 0 0 0 1 0 0]' + [0 0 0 0 0 0 0 0 0 0 -1]';
        
        
        
        [P2m0Ortho] = Orthogonal (P2m0Ket, InitialBra);
        [P2m1Ortho] = Orthogonal (P2m1Ket, InitialBra);
        [P1m2Ortho] = Orthogonal (P1m2Ket, InitialBra);
        [P2m2Ortho] = Orthogonal (P2m2Ket, InitialBra);
        [P1m3Ortho] = Orthogonal (P1m3Ket, InitialBra);
        [P2m3Ortho] = Orthogonal (P2m3Ket, InitialBra);
        [P1m4Ortho] = Orthogonal (P1m4Ket, InitialBra);
        [P2m4Ortho] = Orthogonal (P2m4Ket, InitialBra);
        [P1m5Ortho] = Orthogonal (P1m5Ket, InitialBra);
        [P2m5Ortho] = Orthogonal (P2m5Ket, InitialBra);
        [P1m6Ortho] = Orthogonal (P1m6Ket, InitialBra);
        [P2m6Ortho] = Orthogonal (P2m6Ket, InitialBra);
        [P1m7Ortho] = Orthogonal (P1m7Ket, InitialBra);
        [P2m7Ortho] = Orthogonal (P2m7Ket, InitialBra);
        [P1m8Ortho] = Orthogonal (P1m8Ket, InitialBra);
        [P2m8Ortho] = Orthogonal (P2m8Ket, InitialBra);
        
        if P2m0Ortho ==1
            VFinal(bra,ket) = VFinal(bra,ket) + P2m0Const;
        end
        if P2m1Ortho ==1
            VFinal(bra,ket) = VFinal(bra,ket) + P2m1Const;
        end
        if P1m2Ortho ==1
            VFinal(bra,ket) = VFinal(bra,ket) + P1m2Const;
        end
        if P2m2Ortho ==1
            VFinal(bra,ket) = VFinal(bra,ket) + P2m2Const;
        end
        if P1m3Ortho ==1
            VFinal(bra,ket) = VFinal(bra,ket) + P1m3Const;
        end
        if P2m3Ortho ==1
            VFinal(bra,ket) = VFinal(bra,ket) + P2m3Const;
        end
        if P1m4Ortho ==1
            VFinal(bra,ket) = VFinal(bra,ket) + P1m4Const;
        end
        if P2m4Ortho ==1
            VFinal(bra,ket) = VFinal(bra,ket) + P2m4Const;
        end
        if P1m5Ortho ==1
            VFinal(bra,ket) = VFinal(bra,ket) + P1m5Const;
        end
        if P2m5Ortho ==1
            VFinal(bra,ket) = VFinal(bra,ket) + P2m5Const;
        end
        if P1m6Ortho ==1
            VFinal(bra,ket) = VFinal(bra,ket) + P1m6Const;
        end
        if P2m6Ortho ==1
            VFinal(bra,ket) = VFinal(bra,ket) + P2m6Const;
        end
        if P1m7Ortho ==1
            VFinal(bra,ket) = VFinal(bra,ket) + P1m7Const;
        end
        if P2m7Ortho ==1
            VFinal(bra,ket) = VFinal(bra,ket) + P2m7Const;
        end
        if P1m8Ortho ==1
            VFinal(bra,ket) = VFinal(bra,ket) + P1m8Const;
        end
    end
end
VFinal