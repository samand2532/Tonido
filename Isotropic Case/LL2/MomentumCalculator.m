%%%%%Mom Calculator
clc; clear all;

KPos = [ 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22; %Matlab position
    0 0 0 0 0 0 0 0 0 0  0  1  1  1  1  1  1  1  1  1  1  1 ; %n value
    -1 0 1 2 3 4 5 6 7 8  9 -1  0  1  2  3  4  5  6  7  8  9]; %Mt value

AllBras = csvread('LuisBasisLL2.csv');

for Count = 1:length(AllBras)
    
    MatA(Count,1) = dot(AllBras(Count,:),KPos(3,:));
end
