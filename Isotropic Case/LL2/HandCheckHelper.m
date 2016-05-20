%%%% HAnd checker
clc; clear all;


KPos = [ 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22; %Matlab position
    0 0 0 0 0 0 0 0 0 0  0  1  1  1  1  1  1  1  1  1  1  1 ; %n value
    -1 0 1 2 3 4 5 6 7 8  9 -1  0  1  2  3  4  5  6  7  8  9]; %Mt value


k1 = 2 ;
k2 = 14;
l1 = 2;
l2 = 14;
mk1 = KPos(3,k1);
mk2 = KPos(3,k2);
ml1 = KPos(3,l1);
ml2 = KPos(3,l2);
nk1 = KPos(2,k1);
nk2 = KPos(2,k2);
nl1 = KPos(2,l1);
nl2 = KPos(2,l2);

absMt = abs(mk1)+abs(mk2)+abs(ml1)+abs(ml2);
Mom = 1/2^(absMt/2)
Root = sqrt(1/(factorial(nk1+abs(mk1))*factorial(nk2+ abs(mk2))*factorial(nl1+ abs(ml1))*factorial(nl2 + abs(ml2))))
I2 = factorial(absMt/2)
