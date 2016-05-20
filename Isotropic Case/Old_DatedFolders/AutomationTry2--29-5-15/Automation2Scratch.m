clear all;

Bra = [ 5 0 0 0 1]'; Ket = [ 5 0 0 0 1]';
g = 1;
m1 = 4; m2 = 0; m3 =4; m4 = 0;
fail = 'fail!';
m1Vect = zeros(5,1); m2Vect = zeros(5,1); m3Vect = zeros(5,1); m4Vect = zeros(5,1);
m1Vect((m1+1),1)=-1;m2Vect((m2+1),1)=-1;m3Vect((m3+1),1)=-1;m4Vect((m4+1),1)=-1;
BigConst = (factorial(m1 + m2)/ sqrt(factorial(m1)*factorial(m2)*factorial(m3)*factorial(m4))) * (1/ 2^(m1+m2+1))*(g/pi);
Orth=0;

Constm4 = sqrt(Ket((m4+1),1));
trans4 = Ket+m4Vect;
Constm3 = sqrt(trans4((m3+1),1));
FinalKet = trans4+m3Vect;

Constm1 = sqrt(Bra((m1+1),1));
trans1 = Bra + m1Vect;
Constm2 = sqrt(trans1((m2+1),1));
FinalBra = trans1 + m2Vect;
    
if (FinalBra(1,1) == FinalKet(1,1)) & (FinalBra(2,1) == FinalKet(2,1)) & (FinalBra(3,1) == FinalKet(3,1)) & (FinalBra(4,1) == FinalKet(4,1));
    Orth = 1;
else Orth = 0;
end



if Orth == 1;
    ConstComb = Constm4*Constm3*Constm2*Constm1;
    BigConst;
    FinalConst = ConstComb * BigConst
else disp(fail)
end











