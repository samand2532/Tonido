%%% pulls out even basis from all kets for LuisV.mat
clc; clear all; tic

KPos = [ 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22; %Matlab position
    0 0 0 0 0 0 0 0 0 0  0  1  1  1  1  1  1  1  1  1  1  1 ; %n value
    -1 0 1 2 3 4 5 6 7 8  9 -1  0  1  2  3  4  5  6  7  8  9]; %Mt value

AllBras = csvread('LuisBasisLL2.csv');

load('LuisV.mat');
VMatnoA322 = struct2array(load('LuisV.mat'));
AllBras = csvread('LuisBasisLL2.csv');
a=0;
for i = 1:length(AllBras)
        
        if mod(dot(AllBras(i,:),KPos(3,:)),2) == 0
            a = a +1;
            MatA(a,:) = VMatnoA322(i,:);
            
            
        end
end
a = 0;
for i = 1:length(AllBras)
    
    if mod(dot(AllBras(i,:),KPos(3,:)),2) == 0
    a = a + 1;
    MatB(:,a) = MatA(:,i);
    end
end








toc