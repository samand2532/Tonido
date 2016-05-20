clc; clear all;
tic
A = 0.03;
g = 1;
N = 6;
LaundauLevel = 2; %%% 1 = LLL, 2 = LL2

KPos = [ 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19; %Matlab position
         0 0 1 0 1 0 1 0 1 0  1  0  1  0  1  0  1  0  1 ; %n value
        -1 0 0 1 1 2 2 3 3 4  4  5  5  6  6  7  7  8  8]; %Mt value

%% Section loads N and L mat, also calls basis vector
if LaundauLevel == 1;
    AllBras = csvread('NewLLL.csv');
    AllKets = AllBras';
    NMat = eye(length(AllBras)).*N;
    LMat = zeros(length(AllBras));
    UMat = zeros(length(AllBras));
    AMat = zeros(length(AllBras));
    for LBra = 1:length(AllBras)
        for LKet = 1:length(AllKets)
            InitialLBra = AllBras(LBra,:);
            InitialLKet = AllKets(:,LKet);
            if InitialLBra == InitialLKet'
                LMat(LBra,LKet) = dot(InitialLBra,KPos(3,:));
            end
        end
    end
elseif LaundauLevel == 2;
    
    AllBras = csvread('NewLL2.csv');
    AllKets = AllBras';
    NMat = eye(length(AllBras)).*N;
    LMat = zeros(length(AllBras));
    UMat = zeros(length(AllBras));
    AMat = zeros(length(AllBras));
    for LBra = 1:length(AllBras)
        for LKet = 1:length(AllKets)
            InitialLBra = AllBras(LBra,:);
            InitialLKet = AllKets(:,LKet);
            if InitialLBra == InitialLKet'
                LMat(LBra,LKet) = dot(InitialLBra,KPos(3,:));
            end
        end
    end
end
%% Interaction Part
OVERLAPTEST = 0;
IntDelta = csvread('NEWINTdelta.csv');
l2Vect = zeros(19,1);
l1Vect = zeros(19,1);
k1Vect = zeros(19,1);
k2Vect = zeros(19,1);
for UBra = 1:length(AllBras)
    for UKet = 1:length(AllBras)
        
        InitialBra = AllBras(UBra,:);
        InitialKet = AllKets(:,UKet);
        
        for DeltaScroll = 1:length(IntDelta)
            
            l2Vect = l2Vect.*0;
            l1Vect = l1Vect.*0;
            k1Vect = k1Vect.*0;
            k2Vect = k2Vect.*0;
            
            k1 = IntDelta(DeltaScroll,1);
            k2 = IntDelta(DeltaScroll,2);
            l1 = IntDelta(DeltaScroll,3);
            l2 = IntDelta(DeltaScroll,4);
            
            nk1 = IntDelta(DeltaScroll,5);
            nk2 = IntDelta(DeltaScroll,6);
            nl1 = IntDelta(DeltaScroll,7);
            nl2 = IntDelta(DeltaScroll,8);
            
            mk1 = IntDelta(DeltaScroll,9);
            mk2 = IntDelta(DeltaScroll,10);
            ml1 = IntDelta(DeltaScroll,11);
            ml2 = IntDelta(DeltaScroll,12);
            
            l2Vect(l2,1) = -1;
            l1Vect(l1,1) = -1;
            k2Vect(k2,1) = 1;
            k1Vect(k1,1) = 1;
            
            l2Const = sqrt(InitialKet(l2,1));
            l2trans = InitialKet + l2Vect;
            if any(l2trans < 0)
                %disp('l2<0!')
                continue
            end
            l1Const = sqrt(l2trans(l1,1));
            l1trans = l2trans + l1Vect;
            if any(l1trans < 0)
               % disp('l1<0!')
                continue
            end       
            k2Const = sqrt(l1trans(k2,1)+1);
            k2trans = l1trans + k2Vect;
            k1Const = sqrt(k2trans(k1,1)+1);
            FinalKet = k2trans + k1Vect;
            
            if FinalKet == InitialBra'
                OVERLAPTEST = OVERLAPTEST +1;
                OpConst = l2Const*l1Const*k2Const*k1Const;
                %I2func = factorial((abs(mk1)+abs(mk2)+abs(ml1)+abs(ml2))/2);
                I2func = IntDelta(DeltaScroll,15);
                RootTerm = IntDelta(DeltaScroll,14);
                Momterm = IntDelta(DeltaScroll,13);
                
                UMat(UBra,UKet) = UMat(UBra,UKet) + (RootTerm*I2func*Momterm*OpConst);
                
                %                 input(num2str(DeltaScroll));
                %                 input(num2str(UBra));
                %                 input(num2str(UKet));
%                 OVERLAPTESTmat(OVERLAPTEST,1) = UBra;
%                 OVERLAPTESTmat(OVERLAPTEST,2) = UKet;
%                 OVERLAPTESTmat(OVERLAPTEST,3) = DeltaScroll;
                
            end
        end
    end
end
%% Aniso Part
if A > 0
    AnisoDelta = csvread('NewAnisoDeltaIncConst.csv');
    %AnisoDelta = csvread('AnisoDeltaWithConst.csv');
    k1Vect = zeros(19,1);
    k2Vect = zeros(19,1);
    
    for ABra = 1:length(AllBras)
        for AKet = 1:length(AllBras)
            
            InitialBra = AllBras(ABra,:);
            InitialKet = AllKets(:,AKet);
            
            for ADeltaScroll = 1:length(AnisoDelta)
                
                k1Vect = k1Vect.*0;
                k2Vect = k2Vect.*0;
                
                k1 = AnisoDelta(ADeltaScroll,1);
                k2 = AnisoDelta(ADeltaScroll,2);
                nk1 = AnisoDelta(ADeltaScroll,3);
                nk2 = AnisoDelta(ADeltaScroll,4);
                mk1 = AnisoDelta(ADeltaScroll,5);
                mk2 = AnisoDelta(ADeltaScroll,6);
                
                k1Vect(k1,1) = 1;
                k2Vect(k2,1) = -1;
                
                k2Const = sqrt(InitialKet(k2,1));
                k2trans = InitialKet + k2Vect;
                if any(k2trans <0 )
                    continue
                end
                k1Const = sqrt(k2trans(k1,1)+1);
                FinalKet = k2trans + k1Vect;
                
                if FinalKet == InitialBra'
                    OpConst = k2Const*k1Const;
                    %I1func = factorial((abs(mk1)+abs(mk2)+2)/2);
                    I1func = AnisoDelta(ADeltaScroll,8);
                    RootTerm = AnisoDelta(ADeltaScroll,7);
                    
                    AMat(ABra,AKet) = AMat(ABra,AKet) + (RootTerm*OpConst*I1func);
                end
            end
        end
    end
end
%% Graphing Part
ss=0;
for Omega = 0.7:0.0001:0.85
    %Omega = 0.776
    
    
    
    ss=ss+1;
    
    NAll = NMat;
    UAll = round(UMat,10).*((g/(4*pi)));
    LAll = LMat.*(1-Omega);
    
    AAll = 0.5.*AMat.*A;
    
    Total = NAll+UAll+LAll+round(AAll,10);
    Eig = eig(Total);
    
    OmegaIncreaseMatrix(1,ss) = Omega;
    OmegaEnergyMatrix(:,ss) = Eig;
    OmegaMatRestricted = OmegaEnergyMatrix(1:9,:);
    
    [V D] = eig(Total);
end

plot(OmegaIncreaseMatrix,OmegaEnergyMatrix);
xlabel('Omega')
ylabel('<E>')
if LaundauLevel == 1
 ylim([8.2 8.9])
 xlim([0.7 0.8])
elseif LaundauLevel == 2
    ylim([7.7 8.5])
 xlim([0.75 0.85])
end
    
%hold on
toc