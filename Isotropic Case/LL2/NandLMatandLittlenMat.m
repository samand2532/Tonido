clc; clear all; tic

KPos = [ 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22; %Matlab position
         0 0 0 0 0 0 0 0 0 0  0  1  1  1  1  1  1  1  1  1  1  1 ; %n value
        -1 0 1 2 3 4 5 6 7 8  9 -1  0  1  2  3  4  5  6  7  8  9]; %Mt value



%AllBras = csvread('LuisBasisLL2.csv');
AllBras = csvread('EvenBasisLL2.csv');
AllKets = AllBras';
NMat = eye(length(AllBras)).*6;



absLMat = zeros(length(AllBras));
LMat = zeros(length(AllBras));
for CountA = 1:length(AllBras)
    for CountB = 1:length(AllKets)
        InitialBra = AllBras(CountA,:);
        InitialKet = AllKets(:,CountB);
    if InitialBra == InitialKet'
        absLMat(CountA,CountB) = dot(AllBras(CountA,:),abs(KPos(3,:)));
        LMat(CountA,CountB) = dot(AllBras(CountA,:),KPos(3,:));
        
    end
    end
end

littlenMat = zeros(length(AllBras));
for CountC = 1:length(AllBras)
    for CountD = 1:length(AllKets)
        InitialBra = AllBras(CountC,:);
        InitialKet = AllKets(:,CountD);
        if InitialBra == InitialKet'
            littlenMat(CountC,CountD) = dot(InitialBra(1,:),KPos(2,:));
        end
    end
end



csvwrite('NMat.csv',NMat);
csvwrite('LMat.csv',LMat);
csvwrite('littlenMat.csv',littlenMat);
csvwrite('absLMat.csv',absLMat);







toc
