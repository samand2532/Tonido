%%%Confluent-hypergeometric function
tic
%kummerU(a,b,z)
for z = -5:1:2000
a = kummerU(0.1,0.2,z);

Tab(z+6,1) = a;


end


plot(real(Tab(:,1)))
toc