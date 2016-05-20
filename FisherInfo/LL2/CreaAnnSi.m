function [ OutVect, Const ] = CreaAnnSi( InVect, P1, P2)
%ann crea op for a^d a, ORDER OF P1...2 is VERY IMPORTANT, DO IT AS
%WOULD ON PAPER!!!!!!

TestVect = InVect;
ann1Vect = zeros(size(TestVect)); 
creVect1 = zeros(size(TestVect)); 
ann1Vect(1,P2) = -1;
creVect1(1,P1) = 1;

aV1 = TestVect + ann1Vect;
if any(aV1 < 0 )
    OutVect = zeros(1,length(TestVect));
    Const = 0;
    return
else ann1Const = sqrt(TestVect(1,P2));
        cV1 = aV1 + creVect1;
        c1Const = sqrt(cV1(1,P1));
    end

OutVect = cV1;
Const = ann1Const*c1Const;

end