function [] = LMatLLL( LLVal )
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here

if LLVal == 1
AllBras = csvread('StatesLLL.csv');
elseif LLVavl ==2
    AllBras = csvread('StatesLLL_and_LL2.csv');
end
AllKets = AllBras';

MomValue = [ -1 0 0 1 1 2 2 3 3 4 4 5 5 6 6 7 7 8 8];

for BraL = 1:length(AllBras)
    for KetL = 1:length(AllKets)
        InitialBraL = AllBras(BraL,:);
        InitialKetL = AllKets(:,KetL);
        
        if InitialBraL == InitialKetL'
            LMat(BraL,KetL) = dot(InitialKetL, MomValue);
        end
        
    end
end
csvwrite('LMatLL2.csv',LMat)

end

