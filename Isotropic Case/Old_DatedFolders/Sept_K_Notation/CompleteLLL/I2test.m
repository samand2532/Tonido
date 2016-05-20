
Mt = 3;
nk1 = 0; nk2 = 0; nl1 =0; nl2 = 0;
mk1= 1; mk2 = 0; ml1 = 2; ml2 = 3;


    fun = @(x) exp(-x).*(x.^(Mt)).* (laguerreL(nk1,abs(mk1),(x/2))) .* (laguerreL(nk2,abs(mk2),(x/2))) .* (laguerreL(nl1,abs(ml1),(x/2))) .* (laguerreL(nl2,abs(ml2),(x/2)));
    
I2func = integral(fun,0,inf)


