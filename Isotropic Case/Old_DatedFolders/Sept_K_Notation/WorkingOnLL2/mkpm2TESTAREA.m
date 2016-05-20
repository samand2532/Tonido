KPos = [ 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19; %Matlab position
          0 0 1 0 1 0 1 0 1 0  1  0  1  0  1  0  1  0  1 ; %n value
         -1 0 0 1 1 2 2 3 3 4  4  5  5  6  6  7  7  8  8]; %Mt value
     
     
     %%%% Looping K to find mk1, mk+/- 2
     aa=0;bb=0;
for kk1 = 1:19
    for kk2 = 1:19
        mk1 = KPos(3,kk1); mk2 = KPos(3,kk2);
        
        if mk1 == (mk2+2)
            aa=aa+1;
            deltaP(aa,1) = kk1;
            deltaP(aa,2) = kk2;
        end
        if mk1 ==(mk2-2)
            bb=bb+1;
            deltaM(bb,1) = kk1;
            deltaM(bb,2) = kk2;
        end
    end
end
delta = [deltaP;deltaM];