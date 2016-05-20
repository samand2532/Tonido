clear all; clc;
tic
AllBras =csvread('StatesLLL.csv');
AllKets = AllBras';

Lmax = 8;
g = 1;
A = 0.03;
N = 6;

KPos = [ 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19; %Matlab position
         0 0 1 0 1 0 1 0 1 0  1  0  1  0  1  0  1  0  1 ; %n value
        -1 0 0 1 1 2 2 3 3 4  4  5  5  6  6  7  7  8  8]; %Mt value
% k(n,m)
%           mt-1   | mt0         | mt1         | mt2         | mt3          | ...
%           (0,-1) | (0,0) (1,0) | (0,1) (1,1) | (0,2) (1,2) |  (0,3) (1,3) | ...
% Mat posit   1        2     3       4     5       6     7        8     9
MomValue = [ -1 0 0 1 1 2 2 3 3 4 4 5 5 6 6 7 7 8 8];
NPosVect = [  0 0 1 0 1 0 1 0 1 0 1 0 1 0 1 0 1 0 1];
%%%%% Plus 2 section

nnP=0;qqP=0;rrP=0;ssP=0;
StorageMatrixAllP = zeros(19^4,4);
for mk1P = 1:19;
    for mk2P = 1:19;
        for ml1P = 1:19;
            for ml2P = 1:19;
                nnP = nnP+1;
                StorageMatrixAllP(nnP,1) = mk1P;
                StorageMatrixAllP(nnP,2) = mk2P;
                StorageMatrixAllP(nnP,3) = ml1P;
                StorageMatrixAllP(nnP,4) = ml2P;
            end
        end
    end
end
mk1ValVectP = zeros(19,1);
mk2ValVectP = zeros(19,1);
ml1ValVectP = zeros(19,1);
ml2ValVectP = zeros(19,1);
for pp = 1:(19^4)
    mk1ValVectP = mk1ValVectP.*0;
    mk2ValVectP = mk2ValVectP.*0;
    ml1ValVectP = ml1ValVectP.*0;
    ml2ValVectP = ml2ValVectP.*0;
    
    mk1ValVectP(StorageMatrixAllP(pp,1),1) = 1;
    mk2ValVectP(StorageMatrixAllP(pp,2),1) = 1;
    ml1ValVectP(StorageMatrixAllP(pp,3),1) = 1;
    ml2ValVectP(StorageMatrixAllP(pp,4),1) = 1;
    
    mk1ValP = dot(mk1ValVectP,MomValue);
    mk2ValP = dot(mk2ValVectP,MomValue);
    ml1ValP = dot(ml1ValVectP,MomValue);
    ml2ValP = dot(ml2ValVectP,MomValue);
    mLP = mk1ValP + mk2ValP;
    mRP = ml1ValP + ml2ValP;
      
    %SMmn = zeros(9040,4);
    if (mLP == (mRP+2) & mLP > 0 ); %& mL <= Lmax)
        qqP = qqP + 1;
        %%%Next 4 lines in terms of Matlab position, equivelent to Mt
        %%%notation but uses vecto position instead.
        SMmnP(qqP,1) = StorageMatrixAllP(pp,1);
        SMmnP(qqP,2) = StorageMatrixAllP(pp,2);
        SMmnP(qqP,3) = StorageMatrixAllP(pp,3);
        SMmnP(qqP,4) = StorageMatrixAllP(pp,4);
    end
end

% following creates n=0 and n=1 tabless
 for rrP = 1:length(SMmnP)
     ssP =ssP+1;
     x1P=0;x2P=0;x3P=0;x4P=0;
     x1P = SMmnP(ssP,1); x2P = SMmnP(ssP,2);x3P=SMmnP(ssP,3);x4P=SMmnP(ssP,4);
     y1P = KPos(2,x1P); y2P = KPos(2,x2P); y3P = KPos(2,x3P); y4P = KPos(2,x4P);
     % next is 0 n value 
     if (y1P + y2P + y3P +y4P) == 0
         APAn0P(ssP,1) = SMmnP(ssP,1);
         APAn0P(ssP,2) = SMmnP(ssP,2);
         APAn0P(ssP,3) = SMmnP(ssP,3);
         APAn0P(ssP,4) = SMmnP(ssP,4);
     end
 end
 APAinc0P = [APAn0P]; zzP=0; ffP=0;zzzP=0;
for kxP = 1:length(APAinc0P)
    zzP=zzP+1;
    if (APAinc0P(zzP,1)) + (APAinc0P(zzP,1)) + (APAinc0P(zzP,1)) + (APAinc0P(zzP,1)) > 0
        APAFinalinc0sP(zzP,1) = APAinc0P(zzP,1);
        APAFinalinc0sP(zzP,2) = APAinc0P(zzP,2);
        APAFinalinc0sP(zzP,3) = APAinc0P(zzP,3);
        APAFinalinc0sP(zzP,4) = APAinc0P(zzP,4);
    end
end
      for fffP = 1:length(APAFinalinc0sP)
        zzzP=zzzP+1;
    if (APAFinalinc0sP(zzzP,1) + APAFinalinc0sP(zzzP,2) + APAFinalinc0sP(zzzP,3) + APAFinalinc0sP(zzzP,4)) >0
        ffP = ffP+1;
        APAno0sP(ffP,1) = APAFinalinc0sP(zzzP,1);
        APAno0sP(ffP,2) = APAFinalinc0sP(zzzP,2);
        APAno0sP(ffP,3) = APAFinalinc0sP(zzzP,3);
        APAno0sP(ffP,4) = APAFinalinc0sP(zzzP,4);
    end
    end
  APAFinalnotP = [APAno0sP];
  xxxP=0;
  for yyyP = 1:length(APAFinalnotP)
      if APAFinalnotP(yyyP,1) == 1 || APAFinalnotP(yyyP,2) == 1 || APAFinalnotP(yyyP,3) == 1 || APAFinalnotP(yyyP,4) == 1;
      else
          xxxP = xxxP+1;
          APAFinalP(xxxP,1) = APAFinalnotP(yyyP,1);
          APAFinalP(xxxP,2) = APAFinalnotP(yyyP,2);
          APAFinalP(xxxP,3) = APAFinalnotP(yyyP,3);
          APAFinalP(xxxP,4) = APAFinalnotP(yyyP,4);
      end
  end
  %%%%%%%%% start of -2 section
  nnM=0;qqM=0;rrM=0;ssM=0;
StorageMatrixAllM = zeros(19^4,4);
for mk1M = 1:19;
    for mk2M = 1:19;
        for ml1M = 1:19;
            for ml2M = 1:19;
                nnM = nnM+1;
                StorageMatrixAllM(nnM,1) = mk1M;
                StorageMatrixAllM(nnM,2) = mk2M;
                StorageMatrixAllM(nnM,3) = ml1M;
                StorageMatrixAllM(nnM,4) = ml2M;
            end
        end
    end
end
mk1ValVectM = zeros(19,1);
mk2ValVectM = zeros(19,1);
ml1ValVectM = zeros(19,1);
ml2ValVectM = zeros(19,1);
for pp = 1:(19^4)
    mk1ValVectM = mk1ValVectM.*0;
    mk2ValVectM = mk2ValVectM.*0;
    ml1ValVectM = ml1ValVectM.*0;
    ml2ValVectM = ml2ValVectM.*0;
    
    mk1ValVectM(StorageMatrixAllM(pp,1),1) = 1;
    mk2ValVectM(StorageMatrixAllM(pp,2),1) = 1;
    ml1ValVectM(StorageMatrixAllM(pp,3),1) = 1;
    ml2ValVectM(StorageMatrixAllM(pp,4),1) = 1;
    
    mk1ValM = dot(mk1ValVectM,MomValue);
    mk2ValM = dot(mk2ValVectM,MomValue);
    ml1ValM = dot(ml1ValVectM,MomValue);
    ml2ValM = dot(ml2ValVectM,MomValue);
    mLM = mk1ValM + mk2ValM;
    mRM = ml1ValM + ml2ValM;
      
    %SMmn = zeros(9040,4);
    if (mLM == (mRM-2) & mLM > 0 ); %& mL <= Lmax)
        qqM = qqM + 1;
        %%%Next 4 lines in terms of Matlab position, equivelent to Mt
        %%%notation but uses vecto position instead.
        SMmnM(qqM,1) = StorageMatrixAllM(pp,1);
        SMmnM(qqM,2) = StorageMatrixAllM(pp,2);
        SMmnM(qqM,3) = StorageMatrixAllM(pp,3);
        SMmnM(qqM,4) = StorageMatrixAllM(pp,4);
    end
end

% following creates n=0 and n=1 tabless
 for rrM = 1:length(SMmnM)
     ssM =ssM+1;
     x1M=0;x2M=0;x3M=0;x4M=0;
     x1M = SMmnM(ssM,1); x2M = SMmnM(ssM,2);x3M=SMmnM(ssM,3);x4M=SMmnM(ssM,4);
     y1M = KPos(2,x1M); y2M = KPos(2,x2M); y3M = KPos(2,x3M); y4M = KPos(2,x4M);
     % next is 0 n value 
     if (y1M + y2M + y3M +y4M) == 0
         APAn0M(ssM,1) = SMmnM(ssM,1);
         APAn0M(ssM,2) = SMmnM(ssM,2);
         APAn0M(ssM,3) = SMmnM(ssM,3);
         APAn0M(ssM,4) = SMmnM(ssM,4);
     end
 end
 APAinc0M = [APAn0M]; zzM=0; ffM=0;zzzM=0;
for kxM = 1:length(APAinc0M)
    zzM=zzM+1;
    if (APAinc0M(zzM,1)) + (APAinc0M(zzM,1)) + (APAinc0M(zzM,1)) + (APAinc0M(zzM,1)) > 0
        APAFinalinc0sM(zzM,1) = APAinc0M(zzM,1);
        APAFinalinc0sM(zzM,2) = APAinc0M(zzM,2);
        APAFinalinc0sM(zzM,3) = APAinc0M(zzM,3);
        APAFinalinc0sM(zzM,4) = APAinc0M(zzM,4);
    end
end
      for fffM = 1:length(APAFinalinc0sM)
        zzzM=zzzM+1;
    if (APAFinalinc0sM(zzzM,1) + APAFinalinc0sM(zzzM,2) + APAFinalinc0sM(zzzM,3) + APAFinalinc0sM(zzzM,4)) >0
        ffM = ffM+1;
        APAno0sM(ffM,1) = APAFinalinc0sM(zzzM,1);
        APAno0sM(ffM,2) = APAFinalinc0sM(zzzM,2);
        APAno0sM(ffM,3) = APAFinalinc0sM(zzzM,3);
        APAno0sM(ffM,4) = APAFinalinc0sM(zzzM,4);
    end
    end
  APAFinalnotM = [APAno0sM];
  xxxM=0;
  for yyyM = 1:length(APAFinalnotM)
      if APAFinalnotM(yyyM,1) == 1 || APAFinalnotM(yyyM,2) == 1 || APAFinalnotM(yyyM,3) == 1 || APAFinalnotM(yyyM,4) == 1;
      else
          xxxM = xxxM+1;
          APAFinalM(xxxM,1) = APAFinalnotM(yyyM,1);
          APAFinalM(xxxM,2) = APAFinalnotM(yyyM,2);
          APAFinalM(xxxM,3) = APAFinalnotM(yyyM,3);
          APAFinalM(xxxM,4) = APAFinalnotM(yyyM,4);
      end
  end 
    
  APAFinalPM2 = [ APAFinalP; APAFinalM];
  csvwrite('APAFinalPM2.csv',APAFinalPM2);
toc

