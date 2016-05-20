function [] = DeltaInteractionVariableLL( LLVal, LengthVect )
%The input argument specifies if for LLL or LL2
%For the deltainteractionfunction, LLL = 1, LL2incLLL = 2

MomValue = [ -1 0 0 1 1 2 2 3 3 4 4 5 5 6 6 7 7 8 8];  

KPos = [ 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19; %Matlab position
         0 0 1 0 1 0 1 0 1 0  1  0  1  0  1  0  1  0  1 ; %n value
        -1 0 0 1 1 2 2 3 3 4  4  5  5  6  6  7  7  8  8]; %Mt value

nn=0;qq=0;rr=0;ss=0;
StorageMatrixAll = zeros(LengthVect^4,4);
for mk1 = 1:LengthVect;
    for mk2 = 1:LengthVect;
        for ml1 = 1:LengthVect;
            for ml2 = 1:LengthVect;
                nn = nn+1;
                StorageMatrixAll(nn,1) = mk1;
                StorageMatrixAll(nn,2) = mk2;
                StorageMatrixAll(nn,3) = ml1;
                StorageMatrixAll(nn,4) = ml2;
            end
        end
    end
end
mk1ValVect = zeros(LengthVect,1);
mk2ValVect = zeros(LengthVect,1);
ml1ValVect = zeros(LengthVect,1);
ml2ValVect = zeros(LengthVect,1);
for pp = 1:(LengthVect^4)
    mk1ValVect = mk1ValVect.*0;
    mk2ValVect = mk2ValVect.*0;
    ml1ValVect = ml1ValVect.*0;
    ml2ValVect = ml2ValVect.*0;
    
    mk1ValVect(StorageMatrixAll(pp,1),1) = 1;
    mk2ValVect(StorageMatrixAll(pp,2),1) = 1;
    ml1ValVect(StorageMatrixAll(pp,3),1) = 1;
    ml2ValVect(StorageMatrixAll(pp,4),1) = 1;
    
    mk1Val = dot(mk1ValVect,MomValue);
    mk2Val = dot(mk2ValVect,MomValue);
    ml1Val = dot(ml1ValVect,MomValue);
    ml2Val = dot(ml2ValVect,MomValue);
    mL = mk1Val + mk2Val;
    mR = ml1Val + ml2Val;
      
    %SMmn = zeros(9040,4);
    if (mL == mR & mL > 0 ); %& mL <= Lmax)
        qq = qq + 1;
        %%%Next 4 lines in terms of Matlab position, equivelent to Mt
        %%%notation but uses vecto position instead.
        SMmn(qq,1) = StorageMatrixAll(pp,1);
        SMmn(qq,2) = StorageMatrixAll(pp,2);
        SMmn(qq,3) = StorageMatrixAll(pp,3);
        SMmn(qq,4) = StorageMatrixAll(pp,4);
    end
end

% following creates n=0 and n=1 tabless
 for rr = 1:length(SMmn)
     ss =ss+1;
     x1=0;x2=0;x3=0;x4=0;
     x1 = SMmn(ss,1); x2 = SMmn(ss,2);x3=SMmn(ss,3);x4=SMmn(ss,4);
     y1 = KPos(2,x1); y2 = KPos(2,x2); y3 = KPos(2,x3); y4 = KPos(2,x4);
     % next is 0 n value 
     if (y1 + y2 + y3 +y4) == 0
         APAn0(ss,1) = SMmn(ss,1);
         APAn0(ss,2) = SMmn(ss,2);
         APAn0(ss,3) = SMmn(ss,3);
         APAn0(ss,4) = SMmn(ss,4);
     end
     %%%% n = 1 allowed values
     if (y1 + y2 + y3 +y4) == 1
         APAn1(ss,1) = SMmn(ss,1);
         APAn1(ss,2) = SMmn(ss,2);
         APAn1(ss,3) = SMmn(ss,3);
         APAn1(ss,4) = SMmn(ss,4);
     end
 end
 
 if LLVal == 1 
    APAinc0 = [APAn0]; 
 elseif LLVal == 2
     APAinc0 = [APAn0;APAn0];
 end
 
     
 zz=0; ff=0;zzz=0;
for kx = 1:length(APAinc0)
    zz=zz+1;
    if (APAinc0(zz,1)) + (APAinc0(zz,1)) + (APAinc0(zz,1)) + (APAinc0(zz,1)) > 0
        APAFinalinc0s(zz,1) = APAinc0(zz,1);
        APAFinalinc0s(zz,2) = APAinc0(zz,2);
        APAFinalinc0s(zz,3) = APAinc0(zz,3);
        APAFinalinc0s(zz,4) = APAinc0(zz,4);
    end
end
      for fff = 1:length(APAFinalinc0s)
        zzz=zzz+1;
    if (APAFinalinc0s(zzz,1) + APAFinalinc0s(zzz,2) + APAFinalinc0s(zzz,3) + APAFinalinc0s(zzz,4)) >0
        ff = ff+1;
        APAno0s(ff,1) = APAFinalinc0s(zzz,1);
        APAno0s(ff,2) = APAFinalinc0s(zzz,2);
        APAno0s(ff,3) = APAFinalinc0s(zzz,3);
        APAno0s(ff,4) = APAFinalinc0s(zzz,4);
    end
    end
  APAFinalnot = [ 2 2 2 2; APAno0s];
  xxx=0;
  for yyy = 1:length(APAFinalnot)
      if APAFinalnot(yyy,1) == 1 || APAFinalnot(yyy,2) == 1 || APAFinalnot(yyy,3) == 1 || APAFinalnot(yyy,4) == 1;
      else
          xxx = xxx+1;
          APAFinal(xxx,1) = APAFinalnot(yyy,1);
          APAFinal(xxx,2) = APAFinalnot(yyy,2);
          APAFinal(xxx,3) = APAFinalnot(yyy,3);
          APAFinal(xxx,4) = APAFinalnot(yyy,4);
      end
  end
  
  %%%This section tests if the delta here is the same as the working version
  
  %JTest(:,:) = APAFinal(:,:) -2;
  %JTest1(:,:) = JTest(:,:)/2        
    

csvwrite('DeltaInteractionVariableLL.csv',APAFinal)

end

