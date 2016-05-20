clc; clear all;

AllBras = csvread('NewLL2.csv');
AllKets = AllBras';

TestKet = zeros(1,19);
TestKet(1,1) = 0;
TestKet(1,2) = 0;
TestKet(1,3) = 0;
TestKet(1,4) = 5;
TestKet(1,5) = 0;
TestKet(1,6) = 0;
TestKet(1,7) = 0;
TestKet(1,8) = 0;
TestKet(1,9) = 1;
TestKet(1,10) = 0;
TestKet(1,11) = 0;
TestKet(1,12) = 0;
TestKet(1,13) = 0;
TestKet(1,14) = 0;
TestKet(1,15) = 0;
TestKet(1,16) = 0;
TestKet(1,17) = 0;
TestKet(1,18) = 0;
TestKet(1,19) = 0;
TestKet
for CountA = 1:length(AllKets)
    if TestKet' == AllKets(:,CountA)
        CountA
        break
    end
end
