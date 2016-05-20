function [] = NMat( N, LengthAllBras )
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

NMat = eye(LengthAllBras).*N;
csvwrite('NMatLL2Final.csv',NMat);

end

