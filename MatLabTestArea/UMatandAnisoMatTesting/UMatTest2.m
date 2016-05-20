clc; clear all;
Nope = 0;
KPos = [ 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19; %Matlab position
         0 0 1 0 1 0 1 0 1 0  1  0  1  0  1  0  1  0  1 ; %n value
        -1 0 0 1 1 2 2 3 3 4  4  5  5  6  6  7  7  8  8]; %Mt value
    
AllBras = csvread('NewLL2.csv');
AllKets = AllBras';
UMat = csvread('UMat7-12.csv');

for CountZ = 145%:length(AllBras)
    if CountZ > 1
        input('Next? (Enter)')
     end
    TestKetNum = CountZ;
    TestKet = AllKets(:,TestKetNum);
    TestKetMom = dot(TestKet,KPos(3,:));


KetVect = zeros(length(AllBras),1);
KetVect(TestKetNum,1) = 1;

OutPut = UMat*KetVect;
OutputPosition = OutPut(:,1)>0;

ComparisonMat(1,1) = TestKetNum;
ComparisonMat(2,1) = TestKetMom;

I = find(OutputPosition);
a = 1;
for CountZ1 = 1:length(I)
    a = a +1;
    OutKet = AllKets(:,I(CountZ1));
    OutMom = dot(OutKet,KPos(3,:));
    ComparisonMat(1,a) = I(CountZ1);
    ComparisonMat(2,a) = OutMom;
    if any(ComparisonMat(2,2:end) ~= ComparisonMat(2,1))
        disp('error')
        I(CountZ1)
    end
    ComparisonMat(3,a) = dot(OutKet,KPos(2,:));
    OutPutVal = I(CountZ1);
    ComparisonMat(4,a) = OutPut(OutPutVal,1);
    
    
end
ComparisonMat
clear a; clear ComparisonMat
end
