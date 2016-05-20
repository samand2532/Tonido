%byhandcellbycell2
clear all;
hbar = 6.63e-34/(2*pi);
g=1;
OmegaIncreaseMatrix = zeros(1,250);
OmegaEnergyMatrix = zeros(2,250);
OmegaEnergyMatrixGround = zeros(2,250);
OmegaEnergyMatrixGround2 = zeros(2,250);
OmegaEnergyMatrixGround3 = zeros(2,250);
n=0;
Wp=0.75;
N=6;

for Omega = 0.6:0.001:0.85;
    n=n+1;
    Ground = N + N*(N-1)/(4*pi);
    First = N + (1-Omega)+ N*(N-1)/(4*pi);
    sec = N + (1-Omega)* 2 + N*((N-1)/(4*pi));
    thir = N + (1-Omega)* 3 + N*((N-1)/(4*pi))
    
    OmegaIncreaseMatrix(1,n) = Omega;
    OmegaEnergyMatrix(:,n) = First;
    OmegaEnergyMatrixGround(:,n) = Ground;
    OmegaEnergyMatrix2(:,n) = sec;
    OmegaEnergyMatrix3(:,n) = thir;
    

end

plot(OmegaIncreaseMatrix,OmegaEnergyMatrix, OmegaIncreaseMatrix, OmegaEnergyMatrixGround, OmegaIncreaseMatrix , OmegaEnergyMatrix2 ,OmegaIncreaseMatrix, OmegaEnergyMatrix3 )
xlabel('Omega')
ylabel('<E>')
ylim([8.0 9.2])
%xlim([0.7 0.8])
