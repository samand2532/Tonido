clc; clear all;

LuisMat = csvread('LuisU.csv');
MyMat = csvread('MyUMat.csv');


DifferenceMat = LuisMat == MyMat;
invDifferenceMat = not(DifferenceMat);

for row = 1:322
    for col = 1:322

        if (invDifferenceMat(row,col) == 1);
            MatA(row,1) = row;
            MatA(row,2) = col;
        end
    end
end
a = 0;
for ZerosRm = 1:length(MatA)
    if MatA(ZerosRm,1) == 0
    else a = a +1;
        MatB(a,:) = MatA(ZerosRm,:);
    end
end

