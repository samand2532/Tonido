%byhandcellbycell2
clear all;
hbar = 6.63e-34/(2*pi);
g=1;
OmegaIncreaseMatrix = zeros(1,250);
OmegaEnergyMatrix0 = zeros(2,250);
OmegaEnergyMatrix2 = zeros(2,250);
OmegaEnergyMatrix4 = zeros(5,250);
n=0;
Wp=0.75;

for Omega = 0.6:0.001:0.85;
    n=n+1;
    
    
    L0 = [0];
    N0 = [6];
    U0 = [ (30/2 * g/(2*pi))];
    
    L2 = [2*(1-Omega) 0; 0 2*(1-Omega)];
    N2 = eye(2).* 6;
    U2 = [ 29/2 *(g/(2*pi))   (sqrt(5)/2)*(g/(2*pi));
         (sqrt(5)/2)*(g/(2*pi))  (25/2)*(g/(2*pi)) ];
     
    L4 = eye(5).*(4*(1-Omega));
    N4 = eye(5) .* 6;
    U4 = [ 1.69 0.089 0.0771 0 0; 0.089 1.8303 0.0689 0.194 0; 0.0771 0.0689 1.6512 0.2251 0; 0 0.1949 0.2251 1.989 0.3376; 0 0 0 0.3376 1.9099];
    
    L6 = eye(11).*(6*(1-Omega));
    N6 = eye(11) .* 6;
    U6 = [ 1.6164 0.0272 0.0431 0.0000 0.0352 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000;...
           0.0272 1.6612 0.0472 0.0629 0.0385 0.0629 0.0000 0.0000 0.0000 0.0000 0.0000;...
           0.0431 0.0472 1.4274 0.1592 0.0609 0.0796 0.0000 0.1194 0.0000 0.0000 0.0000;...
           0.0000 0.0629 0.1592 1.6711 0.0000 0.0995 0.1194 0.0000 0.0597 0.0000 0.0000;...
           0.0352 0.0385 0.0609 0.0000 1.3230 0.1949 0.0000 0.0000 0.0000 0.0000 0.0000;...
           0.0000 0.0629 0.0796 0.0995 0.1949 1.6114 0.2387 0.1194 0.2387 0.0000 0.0000;...
           0.0000 0.0000 0.0000 0.1194 0.0000 0.2387 1.6711 0.0000 0.1194 0.1949 0.0000;...
           0.0000 0.0000 0.1194 0.0000 0.0000 0.1194 0.0000 1.3727 0.2387 0.0000 0.0000;...
           0.0000 0.0000 0.0000 0.0597 0.0000 0.2387 0.1194 0.2387 1.7308 0.3898 0.0000;...
           0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.1949 0.0000 0.3898 1.6711 0.3082;...
           0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.3082 1.1937];
   
    Lmat0 = L0 ;
    Nmat0 = N0;
    Umat0 =  U0;
    Totmat0 = Lmat0+Nmat0+Umat0;
    Eig0 = eig(Totmat0);
    
    Lmat2 = L2 ;
    Nmat2 = N2;
    Umat2 =  U2;
    Totmat2 = Lmat2+Nmat2+Umat2;
    Eig2 = eig(Totmat2);
    
    Lmat4 = L4 ;
    Nmat4 = N4;
    Umat4 =  U4;
    Totmat4 = Lmat4+Nmat4+Umat4;
    Eig4 = eig(Totmat4);
    
    Lmat6 = L6 ;
    Nmat6 = N6;
    Umat6 =  U6;
    Totmat6 = Lmat6+Nmat6+Umat6;
    Eig6 = eig(Totmat6);
    
    
    
    OmegaIncreaseMatrix(1,n) = Omega;
    OmegaEnergyMatrix0(:,n) = Eig0;
    OmegaEnergyMatrix2(:,n) = Eig2;
    OmegaEnergyMatrix4(:,n) = Eig4;
    OmegaEnergyMatrix6(:,n) = Eig6;
    

end

plot(OmegaIncreaseMatrix,OmegaEnergyMatrix0,OmegaIncreaseMatrix,OmegaEnergyMatrix2, OmegaIncreaseMatrix,OmegaEnergyMatrix4, OmegaIncreaseMatrix,OmegaEnergyMatrix6)
xlabel('Omega')
ylabel('<E>')
ylim([8 9.2])
%xlim([0.7 0.8])
