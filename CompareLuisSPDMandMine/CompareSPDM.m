clc; clear all;

%%%% Loop that spots differences between mine and Luis' SPDM for set A.

LA08 = importdata('LA0_8.mat');
LA089 = importdata('LA0_89.mat');
LA08936 = importdata('LA0_8936.mat');
MA08 = importdata('MA0_8.mat');
MA089 = importdata('MA0_89.mat');
MA08936 = importdata('MA0_8936.mat');
%%%%%%%%% NAME TEST MATS HERE
MyTest = MA089;
LTest = LA089;
Acc = 4; % Acc = number of dp for rounding
%%%%%%%%%
format long
for CountA = 1:22
    for CountB = 1:22  
        MatDiff(CountA,CountB) = MyTest(CountA,CountB) - LTest(CountA,CountB);
        
%         if round(MyTest(CountA,CountB),Acc) ~= round(LTest(CountA,CountB),Acc)
%             disp(['MyVal','  ', num2str(MyTest(CountA,CountB)),'   ','LVal','   ' num2str(LTest(CountA,CountB))])      
%         end
    end
end