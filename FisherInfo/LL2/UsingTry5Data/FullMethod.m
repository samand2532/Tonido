clc; clear all;

syms phi C60 C42 C24 C06 real
format short

Probs = csvread('ProbVaryOmegas.csv');
%phi = 1 %%% CAN I DO THIS???? PHI SHOULD ALWAYS CANCEL SO JUST SET TI TO ARB NUMBER?????
ExcelRow = 1;
C60 = Probs(2,ExcelRow);
C42 = Probs(3,ExcelRow);
C24 = Probs(4,ExcelRow);
C06 = Probs(5,ExcelRow);

PsiTM(1,1) = C60*exp(i*phi*6);
PsiTM(2,1) = C42*exp(i*phi*4);
PsiTM(3,1) = C24*exp(i*phi*2);
PsiTM(4,1) = C06*exp(i*phi*0);

PsiTMPri(:,1) = diff(PsiTM(:,1),phi);

PsiTMPriConj(:,1) = conj(PsiTMPri(:,1));
TT1Orig = PsiTMPriConj.*PsiTMPri;
Fqpa = sum(TT1Orig)
FqpbP1 = sum(PsiTMPriConj.*PsiTM)
conjFqpbP1 = conj(FqpbP1)
Fqpb = FqpbP1.*conjFqpbP1

Fq = double(4*(Fqpa - Fqpb))

