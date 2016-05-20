% OrthoTestL=6
function [Orthogonal] = Orthogonal (FinalBra, FinalKet)
w=0;
Orth=0;
Orthogonal=0;
for q = 1:length(FinalKet)
    w = w+1;
    if FinalBra(w,1) == FinalKet(w,1);
        Orth = Orth +1;
    end
end
if Orth == length(FinalKet)
    Orthogonal = 1;
end
