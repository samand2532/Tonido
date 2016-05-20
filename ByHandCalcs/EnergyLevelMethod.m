%---------Energy Level values from Rod thesis Pg 55
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
    Ground = N + N*g*(N-1)/(4*pi);
    %First = N + (1-Omega)+ g*N*(N-1)/(4*pi);
    sec = N + ((1-Omega)* 2) + N*g*(((2*N)-2-2)/(8*pi));
    four = N + ((1-Omega)* 4) + (N*g*(((2*N)-4-2)/(8*pi)));
    
    
    OmegaIncreaseMatrix(1,n) = Omega;
    %OmegaEnergyMatrix(:,n) = First;
    OmegaEnergyMatrixGround(:,n) = Ground;
    OmegaEnergyMatrix2(:,n) = sec;
    OmegaEnergyMatrix3(:,n) = four;
    

end

plot( OmegaIncreaseMatrix, OmegaEnergyMatrixGround, OmegaIncreaseMatrix , OmegaEnergyMatrix2 ,OmegaIncreaseMatrix, OmegaEnergyMatrix3 )
xlabel('Omega')
ylabel('<E>')
ylim([8.0 9.2])
%xlim([0.7 0.8])
