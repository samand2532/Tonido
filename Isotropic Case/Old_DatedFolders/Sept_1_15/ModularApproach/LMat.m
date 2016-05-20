% Does not include 1-OMEGA TERM

clear all; clc;



AllBras = csvread('AllBrasLLL.csv');
AllKets = AllBras';
MomValue = [ -1 0 1 2 3 4 5 6 7 8 9 10 0]';

for BraL = 1:length(AllBras)
    for KetL = 1:length(AllBras)
        InitialBraL = AllBras(BraL,:);
        InitialKetL = AllKets(:,KetL);
        if InitialBraL == InitialKetL'
            LMatFinal (BraL,KetL) = dot(InitialKetL,MomValue);
        end
    end
end
csvwrite('LMatFinal.csv',LMatFinal)

