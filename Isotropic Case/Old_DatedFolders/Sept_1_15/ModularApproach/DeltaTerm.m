clc; clear all;

AllBras = csvread('AllBrasLLL.csv');
AllKets = AllBras';
Lmax = 8;

StorageMatrixAll = zeros(4^Lmax,5);
nn=0; pp=0; qq=0;
for m1 = 1:(Lmax+1);
    for m2 = 1:(Lmax+1);
        for m3 = 1:(Lmax+1);
            for m4 = 1:(Lmax+1);
                nn = nn+1;
                StorageMatrixAll(nn,1) = m1-1;
                StorageMatrixAll(nn,2) = m2-1;
                StorageMatrixAll(nn,3) = m3-1;
                StorageMatrixAll(nn,4) = m4-1;
            end
        end
    end
end

for pp = 1:(4^Lmax)
    mL = StorageMatrixAll(pp,1)+StorageMatrixAll(pp,2);
    mR = StorageMatrixAll(pp,3)+StorageMatrixAll(pp,4);
    if (mL == mR & mL > 0 & mL <= Lmax)
        qq = qq + 1;
        SM(qq,:) = StorageMatrixAll(pp,:);
    end
end
APA = [0 0 0 0 0; SM];

csvwrite('APA.csv',APA);