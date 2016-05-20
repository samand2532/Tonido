function [ CreaConst, NewVect ] = CreaOp( Pos, InVect )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
OrigBasis = InVect;
CreaVect = zeros(size(OrigBasis));
CreaVect(1,Pos) = 1;
NewVect = OrigBasis+CreaVect;
CreaConst = sqrt(OrigBasis(1,Pos) + 1);
end