
mk1 = 1;
mk2 = 3;
nk1 = 0;
nk2 = 1;



MtTerm = (abs(mk1)+abs(mk2)+2)/2;
fun = @(x) exp(-x).*(x.^(MtTerm)).* (laguerreL(nk1,abs(mk1),(x))) .*...
    (laguerreL(nk2,abs(mk2),(x)));

I2func = integral(fun,0,inf)