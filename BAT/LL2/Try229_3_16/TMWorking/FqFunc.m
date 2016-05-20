function [ Fq ] = untitled2(C60,C42,C24,C06)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
syms phi real

PsiTM(1,1) = C60*exp(i*phi*6);
PsiTM(2,1) = C42*exp(i*phi*4);
PsiTM(3,1) = C24*exp(i*phi*2);
PsiTM(4,1) = C06*exp(i*phi*0);

PsiTMPri(:,1) = diff(PsiTM(:,1),phi);

PsiTMPriConj(:,1) = conj(PsiTMPri(:,1));
TT1Orig = PsiTMPriConj.*PsiTMPri;
Fqpa = sum(TT1Orig);
FqpbP1 = sum(PsiTMPriConj.*PsiTM);
conjFqpbP1 = conj(FqpbP1);
Fqpb = FqpbP1.*conjFqpbP1;

Fq = 4*(Fqpa - Fqpb);

end

