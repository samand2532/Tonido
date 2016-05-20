%%%%% Selects only even momenta total basis

clc; clear all; tic

KPos = [ 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22; %Matlab position
    0 0 0 0 0 0 0 0 0 0  0  1  1  1  1  1  1  1  1  1  1  1 ; %n value
    -1 0 1 2 3 4 5 6 7 8  9 -1  0  1  2  3  4  5  6  7  8  9]; %Mt value


AllBras = csvread('LuisBasisLL2.csv');
AllKets = AllBras;

aa = 0;

for CountA = 1:length(AllBras)
    if dot(AllBras(CountA,:),KPos(3,:)) == 0
        aa = aa + 1;
        EvenBasis(aa,:) = AllBras(CountA,:);
    elseif dot(AllBras(CountA,:),KPos(3,:)) == 2
        aa = aa + 1;
        EvenBasis(aa,:) = AllBras(CountA,:);
    elseif dot(AllBras(CountA,:),KPos(3,:)) == 4
        aa = aa + 1;
        EvenBasis(aa,:) = AllBras(CountA,:);
    elseif dot(AllBras(CountA,:),KPos(3,:)) == 6
        aa = aa + 1;
        EvenBasis(aa,:) = AllBras(CountA,:);
    elseif dot(AllBras(CountA,:),KPos(3,:)) == 8
        aa = aa + 1;
        EvenBasis(aa,:) = AllBras(CountA,:);
    end
end

aa





toc