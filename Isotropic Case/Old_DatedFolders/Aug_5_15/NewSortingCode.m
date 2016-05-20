clc; clear all;
Lmax=8;
nn=0;qq=0;
for m1 = -1:Lmax
    for m2 = -1:Lmax
        for m3 = -1:Lmax
            for m4 = -1:Lmax
                for n_i = 0:4
                nn=nn+1;
                StorageMatrixAll(nn,1) = m1;
                StorageMatrixAll(nn,2) = m2;
                StorageMatrixAll(nn,3) = m3;
                StorageMatrixAll(nn,4) = m4;
                StorageMatrixAll(nn,5) = n_i;
            end
        end
    end
    end
end

for pp=1:(length(StorageMatrixAll))
    mL = StorageMatrixAll(pp,1)+StorageMatrixAll(pp,2);
    mR = StorageMatrixAll(pp,3)+StorageMatrixAll(pp,4);
    if (mL == mR) & (mL <= Lmax) & (mL ~=0)
        qq=qq+1;
        SortingMat(qq,:) = StorageMatrixAll(pp,:);
    end
end
AllPossibleArrangements = [0 0 0 0 0; SortingMat];
%section following removes all occurances of more than 1 lot of -1
remove = AllPossibleArrangements == -1;
remove = sum(remove,2);
jj=0;
for ii = 1:size(remove)
    if (remove(ii)>1)
    else jj = jj + 1; 
        AllPossibleArrangements1(jj,:) = AllPossibleArrangements(ii,:);
    end
end
jj1 = 0;
for jj=1:length(AllPossibleArrangements1)
    
    if ((AllPossibleArrangements1(jj,1) ==-1) || (AllPossibleArrangements1(jj,2) ==-1) || (AllPossibleArrangements1(jj,3) ==-1) || (AllPossibleArrangements1(jj,4) ==-1)) && ((AllPossibleArrangements1(jj,5) > 0))
    
    else jj1 = jj1+1; 
        AllPossibleArrangementsClean(jj1,:) = AllPossibleArrangements1(jj,:);
    end
end
AllPossibleArrangementsClean
