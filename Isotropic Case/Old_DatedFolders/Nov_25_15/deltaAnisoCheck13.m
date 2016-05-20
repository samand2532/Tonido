clc; clear all;
tic


KPos = [ 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19; %Matlab position
         0 0 1 0 1 0 1 0 1 0  1  0  1  0  1  0  1  0  1 ; %n value
        -1 0 0 1 1 2 2 3 3 4  4  5  5  6  6  7  7  8  8]; %Mt value
    
    
    aa = 0;
    for k1 = 1:19
        for k2 = 1:19
            
                    aa = aa + 1;
                    MatA(aa,1) = k1;
                    MatA(aa,2) = k2;
                    MatA(aa,3) = KPos(2,k1);
                    MatA(aa,4) = KPos(2,k2);
                    MatA(aa,5) = KPos(3,k1);
                    MatA(aa,6) = KPos(3,k2);     
        end
    end
    
    bb = 0;
    for CountA = 1:length(MatA)
        if (MatA(CountA,5)) == ((MatA(CountA,6)) +2);
            bb = bb + 1;
            MatB(bb,:) = MatA(CountA,:);
        elseif (MatA(CountA,5)) == ((MatA(CountA,6)) -2);
            bb = bb + 1;
            MatB(bb,:) = MatA(CountA,:);
        end
    end
    %%% Remove multiple n's
    RemoveN = MatB(:,3:4) == 1;
    sumRemoveN = sum(RemoveN,2);
    cc = 0;
    for CountB = 1:length(sumRemoveN)
        if sumRemoveN(CountB) < 2
            cc = cc + 1;
            MatC(cc,:) = MatB(CountB,:);
        end
    end
    
    %%% stop mt = -1 and n=1 occuring together
    
    RemoveNA = MatC(:,3:4) == 1;
    RemoveMt = MatC(:,5:6) == -1;
    
    TotRem = [RemoveNA RemoveMt];
    sumTotRem = sum(TotRem,2);
    
    dd = 0;
    for CountC = 1:length(TotRem)
        if sumTotRem(CountC) < 2
            dd = dd + 1;
            MatD(dd,:) = MatC(CountC,:);
        end
    end
    
    ee = 0;
    for CountD = 1:length(MatD)
        if (MatD(CountD,5) + MatD(CountD,6)) < 9
            ee = ee + 1;
            MatE(ee,:) = MatD(CountD,:);
        end
    end
    
    
    csvwrite('deltaAniso.csv',MatE);
    
    
    
    toc