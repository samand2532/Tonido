% OrthoTestL=6
function [Orthogonal] = Orthogonal (FinalBra, FinalKet)
w=0;
Orth=0;
Orthogonal=0;
for q = 1:length(FinalKet)
    w = w+1;
    if FinalBra(1,w) == FinalKet(1,w);
        Orth = Orth +1;
    end
end
if Orth ==7
    Orthogonal = 1;
end
