clear all; clc;

tic
N = 6;
g = 1;
A = 0.00;

MomVal = [ -1 0 0 1 1 2 2 3 3 4 4 5 5 6 6 7 7 8 8]';

KPos = [ 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19; %Matlab position
    0 0 1 0 1 0 1 0 1 0  1  0  1  0  1  0  1  0  1 ; %n value
    -1 0 0 1 1 2 2 3 3 4  4  5  5  6  6  7  7  8  8]; %Mt value

AllBras = csvread('NewLL2.csv');
AllKets = AllBras';
%%%%%N and L Mats%%%%%%%%%%%%%%%
NMat = eye(length(AllBras)).*N;
IntMat=zeros(size(NMat));
AnisoMat = zeros(size(NMat));

LMat = zeros(length(AllBras));

for LcountA = 1:length(AllBras)
    for LcountB = 1:length(AllBras)
        
        LTestBra = AllBras(LcountA,:);
        LTestKet = AllKets(:,LcountB);
        
        if LTestBra == LTestKet'
            LMat(LcountA,LcountB) = dot(LTestBra,MomVal);
        end
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%Interaction term%%%%%%%%%%

Intdelta = csvread('deltaInt.csv');


%%%%Section creates all of interaction term, excluding first Const and
%%%%OpValues
for IntCountA = 1:length(Intdelta)
    
    k1 = Intdelta(IntCountA,1);
    k2 = Intdelta(IntCountA,2);
    l1 = Intdelta(IntCountA,3);
    l2 = Intdelta(IntCountA,4);
    mk1 = Intdelta(IntCountA,9);
    mk2 = Intdelta(IntCountA,10);
    ml1 = Intdelta(IntCountA,11);
    ml2 = Intdelta(IntCountA,12);
    nk1 = Intdelta(IntCountA,5);
    nk2 = Intdelta(IntCountA,6);
    nl1 = Intdelta(IntCountA,7);
    nl2 = Intdelta(IntCountA,8);
    % I2func savedin Intdelta(IntCountA,13)
    %%% momTerm is(...,14)
    %%% RootTerm (...,15)
    
    halfMt = (abs(mk1)+abs(mk2)+abs(ml1)+abs(ml2))/2;
    fun = @(x) exp(-x).*(x.^(halfMt)).* (laguerreL(nk1,abs(mk1),(x/2))) .*...
        (laguerreL(nk2,abs(mk2),(x/2))) .* (laguerreL(nl1,abs(ml1),(x/2))) .*...
        (laguerreL(nl2,abs(ml2),(x/2)));
    
    
    I2func = integral(fun,0,inf);
    Intdelta(IntCountA,13) = I2func;
    Intdelta(IntCountA,14) = 1/(2^halfMt);
    
    
    RootA =(factorial(nk1)*factorial(nk2)*factorial(nl1)*factorial(nl2));
    RootB = (factorial(nk1+abs(mk1))*factorial(nk2+abs(mk2))*factorial(nl1+abs(ml1))*...
        factorial(nl2+abs(ml2)));
    RootTerm = sqrt(RootA/RootB);
    Intdelta(IntCountA,15) = RootTerm;
end


%%%%% Op trans and const

k1transket = zeros(19,1);
k2transket = zeros(19,1);
l1transket = zeros(19,1);
l2transket = zeros(19,1);



for IntBra = 1:length(AllBras)
    for IntKet = 1:length(AllKets)
        IntBra
        InitialBra = AllBras(IntBra,:);
        InitialKet = AllKets(:,IntKet);
        
        for IntCountB = 1:length(Intdelta)
            
            k1transket = k1transket.*0;
            k2transket = k2transket.*0;
            l1transket = l1transket.*0;
            l2transket = l2transket.*0;
            
            k1 = Intdelta(IntCountB,1);
            k2 = Intdelta(IntCountB,2);
            l1 = Intdelta(IntCountB,3);
            l2 = Intdelta(IntCountB,4);
            
            k1transket(k1,1) = 1;
            k2transket(k2,1) = 1;
            l1transket(l1,1) = -1;
            l2transket(l2,1) = -1;
            
            l2Const = sqrt(InitialKet(l2,1));
            l2trans = InitialKet+l2transket;
            l1Const = sqrt(l2trans(l1,1));
            l1trans = l2trans+l1transket;
            k2Const = sqrt(l1trans(k2,1)+1);
            k2trans = l1trans + k2transket;
            k1Const = sqrt(k2trans(k1,1)+1);
            FinalKet = k2trans + k1transket;
            
            OpConst = l2Const*l1Const*k2Const*k1Const;
            IntConst = Intdelta(IntCountB,13)*Intdelta(IntCountB,14)*Intdelta(IntCountB,15);
            
            if InitialBra' == FinalKet
                IntMat(IntBra,IntKet) = IntMat(IntBra,IntKet) + (OpConst*IntConst);
            end
            
        end
        
    end
end

if A > 0
    
    Anisodelta = csvread('deltaAniso.csv');
    
    for AnisoCountA = 1:length(Anisodelta)
        
        k1 = Anisodelta(AnisoCountA,1);
        k2 = Anisodelta(AnisoCountA,2);
        l1 = Anisodelta(AnisoCountA,3);
        l2 = Anisodelta(AnisoCountA,4);
        mk1 = Anisodelta(AnisoCountA,9);
        mk2 = Anisodelta(AnisoCountA,10);
        ml1 = Anisodelta(AnisoCountA,11);
        ml2 = Anisodelta(AnisoCountA,12);
        nk1 = Anisodelta(AnisoCountA,5);
        nk2 = Anisodelta(AnisoCountA,6);
        nl1 = Anisodelta(AnisoCountA,7);
        nl2 = Anisodelta(AnisoCountA,8);
        
        %%%% rootterm (position 13)
        anisoRootA = (factorial(nk1)*factorial(nk2));
        anisoRootB = factorial(nk1 + abs(mk1))*factorial(nk2 + abs(mk2));
        RootTermAniso = sqrt(anisoRootA/anisoRootB);
        Anisodelta(AnisoCountA,13) = RootTermAniso;
        
        
        funaniso = @(x) exp(-x).*(x^((abs(mk1)+abs(mk2)+2)/2)).*...
            laguerreL(nk1,abs(mk1),x).*laguerreL(nk2,abs(mk2),x);
        I1func = integral(funaniso,0,inf);
        Anisodelta(AnisoCountA,14) = I1func;
    end
    
    for AniBra = 1:length(AllBras)
        for AniKet = 1:length(AllKets)
            
            InitialBra = AllBras(AniBra,:);
            InitialKet = AllKets(:,AniKet);
            
            for AniCountB = 1:length(Anisodelta)
                
                k1transket = k1transket.*0;
                k2transket = k2transket.*0;
                
                k1 = Anisodelta(AniCountB,1);
                k2 = Anisodelta(AniCountB,2);
                
                k1transket(k1,1) = 1;
                k2transket(k2,1) = -1;
                
                k2Const = sqrt(InitialKet(k2,1));
                k2trans = InitialKet + k2transket;
                k1Const = sqrt(k2trans(k1,1)+1);
                FinalKet = k2trans + k1transket;
                                
                AniOpConst = k2Const*k1Const;
                AniConst = Anisodelta(AnisoCountA,13)*Anisodelta(AnisoCountA,14);
                
                if InitialBra' == FinalKet
                    AnisoMat(AniBra,AniKet) = AnisoMat(AniBra,AniKet) + (AniOpConst*AniConst);
                end
                
            end
            
        end
    end
    
else
end

ss=0;
for Omega = 0.6:0.001:1
    ss=ss+1;
    
    NAll = NMat;
    UAll = IntMat.*(g/(4*pi));
    LAll = LMat.*(1-Omega);
    
    VAll = AnisoMat.*A;
    
    Total = NAll+UAll+LAll+VAll;
    Eig = eig(Total);
    
    OmegaIncreaseMatrix(1,ss) = Omega;
    OmegaEnergyMatrix(:,ss) = Eig;
    OmegaMatRestricted = OmegaEnergyMatrix(1:9,:);
    
end

plot(OmegaIncreaseMatrix,OmegaEnergyMatrix);
xlabel('Omega')
ylabel('<E>')
ylim([8 9.2])
xlim([0.6 0.85])


toc


















