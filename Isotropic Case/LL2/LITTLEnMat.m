%%%%% calculate s first term of 3.43
clc; clear all; tic
%AllBras = csvread('LuisBasisLL2.csv');
AllBras = csvread('EvenBasisLL2.csv');
AllKets = AllBras';

KPos = [ 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22; %Matlab position
         0 0 0 0 0 0 0 0 0 0  0  1  1  1  1  1  1  1  1  1  1  1 ; %n value
        -1 0 1 2 3 4 5 6 7 8  9 -1  0  1  2  3  4  5  6  7  8  9]; %Mt value

nMat = zeros(length(AllBras));

for Bras = 1:length(AllBras)
    for Kets = 1:length(AllKets)
        InitialBra = AllBras(Bras,:);
        InitialKet = AllKets(:,Kets);
        
        if InitialBra == InitialKet'
            nMat(Bras,Kets) = dot(InitialBra(1,:),KPos(2,:));
        end
    end
end



csvwrite('nMat.csv',nMat);

toc