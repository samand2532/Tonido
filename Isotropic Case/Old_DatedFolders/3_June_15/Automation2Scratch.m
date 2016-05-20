clear all;

Bra = [ 3 2 1 0 0]'; Ket = [ 3 2 1 0 0]';
syms g p
m1 = 0; m2 = 0; m3 =0; m4 = 0;
fail = 'YOUR COMPUTER VIOLATES MATHS!';
m1Vect = zeros(5,1); m2Vect = zeros(5,1); m3Vect = zeros(5,1); m4Vect = zeros(5,1);
m1Vect((m1+1),1)=-1;m2Vect((m2+1),1)=-1;m3Vect((m3+1),1)=-1;m4Vect((m4+1),1)=-1;
BigConst = (factorial(m1 + m2)/ sqrt(factorial(m1)*factorial(m2)*factorial(m3)*factorial(m4))) * (1/ 2^(m1+m2+1))*(g/p);
Orth=0;OrthoTest =0;

Constm4 = sqrt(Ket((m4+1),1));
trans4 = Ket+m4Vect;
Constm3 = sqrt(trans4((m3+1),1));
FinalKet = trans4+m3Vect;

Constm1 = sqrt(Bra((m1+1),1));
trans1 = Bra + m1Vect;
Constm2 = sqrt(trans1((m2+1),1));
FinalBra = trans1 + m2Vect;

for m = 1:5
    
    if FinalBra(m,1) == FinalKet(m,1);
        Orth = Orth+1;
        if Orth == 5;
            OrthoTest =1;
            if OrthoTest == 1;
                OrthoTest
            else
                disp(fail)
            end
            
        end
    end
end
if OrthoTest == 1;
    ConstComb = Constm4*Constm3*Constm2*Constm1
    BigConst
else disp(fail)
end











