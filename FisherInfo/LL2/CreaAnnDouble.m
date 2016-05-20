function [ OutVect, Const ] = CreaAnnDouble( InVect, P1, P2, P3, P4 )
%ann crea of for a^d a^d a a, ORDER OF P1...4 is VERY IMPORTANT, DO IT AS
%WOULD ON PAPER!!!!!!

TestVect = InVect;
ann1Vect = zeros(size(TestVect)); ann2Vect = zeros(size(TestVect));
creVect1 = zeros(size(TestVect)); creVect2 = zeros(size(TestVect));
ann1Vect(1,P4) = -1; ann2Vect(1,P3) = -1;
creVect1(1,P2) = 1; creVect2(1,P1) = 1;

aV1 = TestVect + ann1Vect;
if any(aV1 < 0 )
    OutVect = zeros(1,length(TestVect));
    Const = 0;
    return
else aV2 = aV1 + ann2Vect;
    ann1Const = sqrt(TestVect(1,P4));
    if any(aV2 < 0)
        OutVect = zeros(1,length(TestVect));
        Const = 0;
        return
    else ann2Const = sqrt(aV2(1,P3)+1);
        cV1 = aV2 + creVect1;
        c1Const = sqrt(cV1(1,P2));
        cV2 = cV1 + creVect2;
        c2Const = sqrt(cV2(1,P1));
    end
end
OutVect = cV2;
Const = c2Const*c1Const*ann2Const*ann1Const;

end

