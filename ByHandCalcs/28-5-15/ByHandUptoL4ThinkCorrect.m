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
    U4 = (g/(pi)).*[ 85/16 (sqrt(5)/8) 0.5*(sqrt(15)/8) 0 0;...
                     (sqrt(5)/8) 0.25*(23) 0.5*(sqrt(3)/4) 0.5*(sqrt(6)/2) 0;...
                     0.5*(sqrt(15)/8) 0.5*(sqrt(3)/4) 0.5*(83/8) 0.5*sqrt(2) 0;...
                     0 0.5*(sqrt(6)/2) 0.5*sqrt(2) 25/4 0.5*((3*sqrt(2))/2);...
                     0 0 0 0.5*((3*sqrt(2))/2) 0.5*12];
    
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
    
    
    
    OmegaIncreaseMatrix(1,n) = Omega;
    OmegaEnergyMatrix0(:,n) = Eig0;
    OmegaEnergyMatrix2(:,n) = Eig2;
    OmegaEnergyMatrix4(:,n) = Eig4;
    

end

plot(OmegaIncreaseMatrix,OmegaEnergyMatrix0,OmegaIncreaseMatrix,OmegaEnergyMatrix2, OmegaIncreaseMatrix,OmegaEnergyMatrix4)
xlabel('Omega')
ylabel('<E>')
ylim([8 9.2])
xlim([0.6 0.85])
