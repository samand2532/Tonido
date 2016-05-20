clc; clear all;
tic


KPos = [ 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19; %Matlab position
         0 0 1 0 1 0 1 0 1 0  1  0  1  0  1  0  1  0  1 ; %n value
        -1 0 0 1 1 2 2 3 3 4  4  5  5  6  6  7  7  8  8]; %Mt value
    
    
    
    
    aa = 0;
    
    for k1 = 1:19
        for k2 = 1:19
            for k3 = 1:19
                for k4 = 1:19
                    aa = aa + 1;
                    MatA(aa,1) = k1;
                    MatA(aa,2) = k2;
                    MatA(aa,3) = k3;
                    MatA(aa,4) = k4;
                    
                    MatA(aa,5) = KPos(2,k1);
                    MatA(aa,6) = KPos(2,k2);
                    MatA(aa,7) = KPos(2,k3);
                    MatA(aa,8) = KPos(2,k4);
                    
                    MatA(aa,9) = KPos(3,k1);
                    MatA(aa,10) = KPos(3,k2);
                    MatA(aa,11) = KPos(3,k3);
                    MatA(aa,12) = KPos(3,k4);
                    
                end
            end
        end
    end
    
    bb = 0;
    for CountB = 1:length(MatA)
        if (MatA(CountB,9) + MatA(CountB,10)) < 9
            bb=bb+1;
            MatB(bb,:) = MatA(CountB,:);
        end
    end
    
    cc = 0;
    for CountC = 1:length(MatB)
        if (MatB(CountC,9)+MatB(CountC,10)) == (MatB(CountC,11)+MatB(CountC,12))
            cc=cc+1;
            MatC(cc,:) = MatB(CountC,:);
        end
    end
    
    %%%% Remove N's
    
    RemoveN = MatC(:,5:8) == 1;
    sumRemoveN = sum(RemoveN,2);
    
    dd = 0;
    for CountD = 1:length(MatC)
        if sumRemoveN(CountD) < 2
            dd = dd+1;
            MatD(dd,:) = MatC(CountD,:);
        end
    end
    
    %%%%% Remove Mt's
    
    RemoveMt = MatD(:,9:12) == -1;
    sumRemoveMt = sum(RemoveMt,2);
    
    ee = 0;
    for CountE = 1:length(MatD)
        if sumRemoveMt(CountE) < 2
            ee = ee +1;
            MatE(ee,:) = MatD(CountE,:);
        end
    end
    
    %%% stop mt and n occuring together
    
    RemoveNa = MatE(:,5:8) == 1;
    RemoveMtA = MatE(:,9:12) == -1
    RemoveAll = [RemoveNa RemoveMtA];
    sumRemAll = sum(RemoveAll,2);
    
    ff = 0;
    
    for CountF = 1:length(MatE)
        if sumRemAll(CountF) < 2
            ff = ff +1;
            MatF(ff,:) = MatE(CountF,:);
        end
    end
    
       
    
    
    toc
    
    
