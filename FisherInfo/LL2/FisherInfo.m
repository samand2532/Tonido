clc; clear all; close all;
tic

KPos = [ 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22; %Matlab position
    0 0 0 0 0 0 0 0 0 0  0  1  1  1  1  1  1  1  1  1  1  1 ; %n value
    -1 0 1 2 3 4 5 6 7 8  9 -1  0  1  2  3  4  5  6  7  8  9]; %Mt value

A = 0.03;
g = 1;

Nt = 6;
%%%% Omega = 0.9375722836016 %best so far for g = 0.4
%%%% Omega = somewhere between 0.82856708053814 for g = 1
%%%% Omega = 0.913847869EnergyLevel
%%%% Omega = 0.8936197 for g = 0.6
%%%% Omega = 0.86697961104 for g = 0.75
%%%% g = 0.2 Cannot find  an omega that causes Bat state, Omega = 0.99999
%%%% isnt high enough.
%%%% Omega = 0.858808166          g = 0.8
%%%% Omega = 0.84326424262349    g = 0.9
OmegaC = 0;

for Omega = 0.8085670805:0.001:0.8485670805
    Omega
    OmegaC = OmegaC + 1;
    Basis = csvread('EvenBasisLL2.csv');
    AllBras = Basis;
    AllKets = AllBras';
    NMat = csvread('NMat.csv');
    absLMat = csvread('absLMat.csv');
    littlenMat = csvread('littlenMat.csv');
    UMat = csvread('MyUMat.csv');
    LMat = csvread('LMat.csv');
    VMatnoA = csvread('VMatnoA.csv');
    
    NAll = NMat;
    LAll = LMat;
    absL = absLMat;
    
    UAll = g.*UMat;
    VAll = VMatnoA.*A;
    
    littlenAll = 2.*littlenMat;
    Total = littlenAll+(absLMat-(Omega*LAll))+(NAll)+UAll+VAll;
    Eig = eig(round(Total,15));
    
    [V,D] = eig(Total);
    
    Dordered = diag(D);
    
    [HamMinVal,HamMinPos] = min(Dordered(:));
    
    [x,y] = size(Basis);
    SPDMmat = zeros(y);
    creaVect = zeros(y,1);
    annVect = zeros(y,1);
    
    OrigCoeff  = V(:,HamMinPos);
    Wavefunc = [OrigCoeff Basis];
    for Kets = 1:x
        for k = 1:y
            for l = 1:y
                C1 =OrigCoeff(Kets,1);
                
                InitialKet = AllKets(:,Kets);
                creaVect = creaVect.*0;
                annVect = annVect.*0;
                creaVect(k,1) = 1;
                annVect(l,1) = -1;
                
                annConst = sqrt(InitialKet(l,1));
                annTrans = InitialKet + annVect;
                creaConst = sqrt(annTrans(k,1)+1);
                FinalKet = annTrans + creaVect;
                
                for Bras = 1:x
                    Bra = AllBras(Bras,:);
                    C2 =OrigCoeff(Bras,1);
                    if Bra == FinalKet'
                        SPDMmat(k,l) = SPDMmat(k,l) + (creaConst*annConst*C1*C2);
                    end
                end
            end
        end
    end
    
    [eigVectSPDM,eigValSPDM] = eig(SPDMmat);
    [MMVal,MMPos] = max(diag(eigValSPDM));
    Vect = eigVectSPDM(:,MMPos);
    [VectMaxVal,VectMaxPos] = max(abs(Vect));
    Vect(VectMaxPos,1) = 0;
    [VectSecMaxVal,VectSecMaxPos] = max(abs(Vect));
    disp(['Most is ' num2str(VectMaxVal) ' at pos A(dag)_' num2str(VectMaxPos)])
    disp(['Sec_Most is ' num2str(VectSecMaxVal) ' at pos A(dag)_' num2str(VectSecMaxPos)])
    
    CreaMat = zeros(2);
    CreaMat(1,1) = VectMaxVal;
    CreaMat(1,2) = VectMaxPos;
    CreaMat(2,1) = VectSecMaxVal;
    CreaMat(2,2) = VectSecMaxPos;
    C2 = CreaMat(1,1);
    C4 = CreaMat(2,1);
    a2 = CreaMat(1,2);
    a4 = CreaMat(2,2);
    CreaMat
    
    %% Second Try NEW FQ METHOD
    %%% FqA
    FqAMat = zeros(192);
    for FqAloop1 = 1:length(Basis)
        InitialVect = Basis(FqAloop1,:);
        [FqAP1Vect, FqAP1Val] = CreaAnnDouble(InitialVect,a2,a2,a2,a2);
        [FqAP2Vect, FqAP2Val] = CreaAnnSi(InitialVect,a2,a2);
        [FqAP3Vect, FqAP3Val] = CreaAnnDouble(InitialVect,a4,a4,a4,a4);
        [FqAP4Vect, FqAP4Val] = CreaAnnSi(InitialVect,a4,a4);
        [FqAP5Vect, FqAP5Val] = CreaAnnDouble(InitialVect,a2,a2,a2,a4);
        [FqAP6Vect, FqAP6Val] = CreaAnnSi(InitialVect,a2,a4);
        [FqAP7Vect, FqAP7Val] = CreaAnnDouble(InitialVect,a2,a4,a2,a4);
        [FqAP8Vect, FqAP8Val] = CreaAnnDouble(InitialVect,a2,a2,a4,a2);
        [FqAP9Vect, FqAP9Val] = CreaAnnDouble(InitialVect,a4,a2,a2,a2);
        [FqAP10Vect, FqAP10Val] = CreaAnnSi(InitialVect,a4,a2);
        [FqAP11Vect, FqAP11Val] = CreaAnnDouble(InitialVect,a2,a4,a2,a4);
        [FqAP12Vect, FqAP12Val] = CreaAnnDouble(InitialVect,a2,a2,a4,a4);
        [FqAP13Vect, FqAP13Val] = CreaAnnDouble(InitialVect,a2,a4,a4,a2);
        [FqAP14Vect, FqAP14Val] = CreaAnnSi(InitialVect,a2,a2);
        [FqAP15Vect, FqAP15Val] = CreaAnnDouble(InitialVect,a4,a2,a2,a4);
        [FqAP16Vect, FqAP16Val] = CreaAnnSi(InitialVect,a4,a4);
        [FqAP17Vect, FqAP17Val] = CreaAnnDouble(InitialVect,a4,a4,a2,a2);
        [FqAP18Vect, FqAP18Val] = CreaAnnDouble(InitialVect,a4,a2,a4,a2);
        [FqAP19Vect, FqAP19Val] = CreaAnnDouble(InitialVect,a2,a4,a4,a4);
        [FqAP20Vect, FqAP20Val] = CreaAnnSi(InitialVect,a2,a4);
        [FqAP21Vect, FqAP21Val] = CreaAnnDouble(InitialVect,a4,a4,a2,a4);
        [FqAP22Vect, FqAP22Val] = CreaAnnDouble(InitialVect,a4,a2,a4,a4);
        [FqAP23Vect, FqAP23Val] = CreaAnnDouble(InitialVect,a4,a4,a4,a2);
        [FqAP24Vect, FqAP24Val] = CreaAnnSi(InitialVect,a4,a2);
        for FqAloop2 = 1:length(Basis)
            if Basis(FqAloop2,:) == FqAP1Vect
                FqAMat(FqAloop2,FqAloop1) = FqAMat(FqAloop2,FqAloop1) + (FqAP1Val * (C2^4)*Wavefunc(FqAloop1,1)*Wavefunc(FqAloop2,1));
            end
            if Basis(FqAloop2,:) == FqAP1Vect
                FqAMat(FqAloop2,FqAloop1) = FqAMat(FqAloop2,FqAloop1) + (FqAP2Val * (C2^4)*Wavefunc(FqAloop1,1)*Wavefunc(FqAloop2,1));
            end
                if Basis(FqAloop2,:) == FqAP3Vect
                FqAMat(FqAloop2,FqAloop1) = FqAMat(FqAloop2,FqAloop1) + (FqAP3Val * (C4^4)*Wavefunc(FqAloop1,1)*Wavefunc(FqAloop2,1));
                end
                if Basis(FqAloop2,:) == FqAP4Vect
                FqAMat(FqAloop2,FqAloop1) = FqAMat(FqAloop2,FqAloop1) + (FqAP4Val * (C4^4)*Wavefunc(FqAloop1,1)*Wavefunc(FqAloop2,1));
                end
                if Basis(FqAloop2,:) == FqAP5Vect
                FqAMat(FqAloop2,FqAloop1) = FqAMat(FqAloop2,FqAloop1) + (FqAP5Val * (C2^3 * C4)*Wavefunc(FqAloop1,1)*Wavefunc(FqAloop2,1));
                end
                if Basis(FqAloop2,:) == FqAP6Vect
                FqAMat(FqAloop2,FqAloop1) = FqAMat(FqAloop2,FqAloop1) + (FqAP6Val * (C2^3 * C4)*Wavefunc(FqAloop1,1)*Wavefunc(FqAloop2,1));
                end
                if Basis(FqAloop2,:) == FqAP7Vect
                FqAMat(FqAloop2,FqAloop1) = FqAMat(FqAloop2,FqAloop1) + (FqAP7Val * (C2^3 * C4)*Wavefunc(FqAloop1,1)*Wavefunc(FqAloop2,1));
                end
                if Basis(FqAloop2,:) == FqAP8Vect
                FqAMat(FqAloop2,FqAloop1) = FqAMat(FqAloop2,FqAloop1) + (FqAP8Val * (C2^3 * C4)*Wavefunc(FqAloop1,1)*Wavefunc(FqAloop2,1));
                end
                if Basis(FqAloop2,:) == FqAP9Vect
                FqAMat(FqAloop2,FqAloop1) = FqAMat(FqAloop2,FqAloop1) + (FqAP9Val * (C2^3 * C4)*Wavefunc(FqAloop1,1)*Wavefunc(FqAloop2,1));
                end
                if Basis(FqAloop2,:) == FqAP10Vect
                FqAMat(FqAloop2,FqAloop1) = FqAMat(FqAloop2,FqAloop1) + (FqAP10Val * (C2^3 * C4)*Wavefunc(FqAloop1,1)*Wavefunc(FqAloop2,1));
                end
                if Basis(FqAloop2,:) == FqAP11Vect
                FqAMat(FqAloop2,FqAloop1) = FqAMat(FqAloop2,FqAloop1) + (FqAP11Val * (C2^2 * C4^2)*Wavefunc(FqAloop1,1)*Wavefunc(FqAloop2,1));
                end
                if Basis(FqAloop2,:) == FqAP12Vect
                FqAMat(FqAloop2,FqAloop1) = FqAMat(FqAloop2,FqAloop1) + (FqAP12Val * (C2^2 * C4^2)*Wavefunc(FqAloop1,1)*Wavefunc(FqAloop2,1));
                end
                if Basis(FqAloop2,:) == FqAP13Vect
                FqAMat(FqAloop2,FqAloop1) = FqAMat(FqAloop2,FqAloop1) + (FqAP13Val * (C2^2 * C4^2)*Wavefunc(FqAloop1,1)*Wavefunc(FqAloop2,1));
                end
                if Basis(FqAloop2,:) == FqAP14Vect
                FqAMat(FqAloop2,FqAloop1) = FqAMat(FqAloop2,FqAloop1) + (FqAP14Val * (C2^2 * C4^2)*Wavefunc(FqAloop1,1)*Wavefunc(FqAloop2,1));
                end
                if Basis(FqAloop2,:) == FqAP15Vect
                FqAMat(FqAloop2,FqAloop1) = FqAMat(FqAloop2,FqAloop1) + (FqAP15Val * (C2^2 * C4^2)*Wavefunc(FqAloop1,1)*Wavefunc(FqAloop2,1));
                end
                if Basis(FqAloop2,:) == FqAP16Vect
                FqAMat(FqAloop2,FqAloop1) = FqAMat(FqAloop2,FqAloop1) + (FqAP16Val * (C2^2 * C4^2)*Wavefunc(FqAloop1,1)*Wavefunc(FqAloop2,1));
                end
                if Basis(FqAloop2,:) == FqAP17Vect
                FqAMat(FqAloop2,FqAloop1) = FqAMat(FqAloop2,FqAloop1) + (FqAP17Val * (C2^2 * C4^2)*Wavefunc(FqAloop1,1)*Wavefunc(FqAloop2,1));
                end
                if Basis(FqAloop2,:) == FqAP18Vect
                FqAMat(FqAloop2,FqAloop1) = FqAMat(FqAloop2,FqAloop1) + (FqAP18Val * (C2^2 * C4^2)*Wavefunc(FqAloop1,1)*Wavefunc(FqAloop2,1));
                end
                if Basis(FqAloop2,:) == FqAP19Vect
                FqAMat(FqAloop2,FqAloop1) = FqAMat(FqAloop2,FqAloop1) + (FqAP19Val * (C2 * C4^3)*Wavefunc(FqAloop1,1)*Wavefunc(FqAloop2,1));
                end
                if Basis(FqAloop2,:) == FqAP20Vect
                FqAMat(FqAloop2,FqAloop1) = FqAMat(FqAloop2,FqAloop1) + (FqAP20Val * (C2 * C4^3)*Wavefunc(FqAloop1,1)*Wavefunc(FqAloop2,1));
                end
                if Basis(FqAloop2,:) == FqAP21Vect
                FqAMat(FqAloop2,FqAloop1) = FqAMat(FqAloop2,FqAloop1) + (FqAP21Val * (C2 * C4^3)*Wavefunc(FqAloop1,1)*Wavefunc(FqAloop2,1));
                end
                if Basis(FqAloop2,:) == FqAP22Vect
                FqAMat(FqAloop2,FqAloop1) = FqAMat(FqAloop2,FqAloop1) + (FqAP22Val * (C2 * C4^3)*Wavefunc(FqAloop1,1)*Wavefunc(FqAloop2,1));
                end
                if Basis(FqAloop2,:) == FqAP23Vect
                FqAMat(FqAloop2,FqAloop1) = FqAMat(FqAloop2,FqAloop1) + (FqAP23Val * (C2 * C4^3)*Wavefunc(FqAloop1,1)*Wavefunc(FqAloop2,1));
                end
                if Basis(FqAloop2,:) == FqAP24Vect
                FqAMat(FqAloop2,FqAloop1) = FqAMat(FqAloop2,FqAloop1) + (FqAP24Val * (C2 * C4^3)*Wavefunc(FqAloop1,1)*Wavefunc(FqAloop2,1));
                
            end
        end
        
    
    end
    %%% FqB
    FqBMat = zeros(192);
    for FqBloop1 = 1:length(Basis)
        InitialVect = Basis(FqBloop1,:);      
        [FqBP1Vect, FqBP1Val] = CreaAnnSi(InitialVect,a2,a2);
        [FqBP2Vect, FqBP2Val] = CreaAnnSi(InitialVect,a2,a4);
        [FqBP3Vect, FqBP3Val] = CreaAnnSi(InitialVect,a4,a2);
        [FqBP4Vect, FqBP4Val] = CreaAnnSi(InitialVect,a4,a4);
        for FqBloop2 = 1:length(Basis)
            if Basis(FqBloop2,:) == FqBP1Vect
                FqBMat(FqBloop2,FqBloop1) = FqBMat(FqBloop2,FqBloop1) + (FqBP1Val * (C2^2)*Wavefunc(FqBloop1,1)*Wavefunc(FqBloop2,1));
            elseif Basis(FqBloop2,:) == FqBP2Vect
                FqBMat(FqBloop2,FqBloop1) = FqBMat(FqBloop2,FqBloop1) + (FqBP2Val * (C2*C4)*Wavefunc(FqBloop1,1)*Wavefunc(FqBloop2,1));
            elseif Basis(FqBloop2,:) == FqBP3Vect
                FqBMat(FqBloop2,FqBloop1) = FqBMat(FqBloop2,FqBloop1) + (FqBP3Val * (C2*C4)*Wavefunc(FqBloop1,1)*Wavefunc(FqBloop2,1));
            elseif Basis(FqBloop2,:) == FqBP4Vect
                FqBMat(FqBloop2,FqBloop1) = FqBMat(FqBloop2,FqBloop1) + (FqBP4Val * (C4*C4)*Wavefunc(FqBloop1,1)*Wavefunc(FqBloop2,1));
            end
        end    
    end

FqBTot = sum(sum(FqBMat.^2))
FqATot = sum(sum(FqAMat))

Fq = 4*(FqATot-FqBTot)
MatTemp(1,OmegaC) = Omega;
MatTemp(2,OmegaC) = Fq;
end
%%%% First Try Fq Method
%% FqA
% FqAMat1 = zeros(192);
% FqAMat2 = zeros(192);
% FqMatAP1 = zeros(192);
% for FqA1loop1 = 1:length(Basis)
% 
%     %%% P1 is for when even number of raising and lowering ON SAME
%     %%% MODE!!!!! ie a^dag _ 2 a_2, (a^dag_2)(a^dag_2)(a_2)(a_2)
% 
%     FqAP1 = (C2^4 * (Basis(FqA1loop1,a2)^2 + (Basis(FqA1loop1,a2)))) + ...
%         (C4^4 * (Basis(FqA1loop1,a4)^2 + (Basis(FqA1loop1,a4)))) +...
%         (C2^2 * C4^2 * (Basis(FqA1loop1,a2)*Basis(FqA1loop1,a4)) + (Basis(FqA1loop1,a2)*Basis(FqA1loop1,a4)) +...
%         (Basis(FqA1loop1,a2)*Basis(FqA1loop1,a4)) + (Basis(FqA1loop1,a2)) );
% 
%     for FqA1loop2 = 1:length(Basis)
%         if Basis(FqA1loop1,:) == Basis(FqA1loop2,:)
%             FqAMat1(FqA1loop1,FqA1loop2) = FqAMat1(FqA1loop1,FqA1loop2) + (FqAP1 * Wavefunc(FqA1loop1,1) * Wavefunc(FqA1loop2,1));
%         end
%     end
% 
%    %%% for terms depending on C^2 _ 2 C^2 _ 4
%         TestVect = Basis(FqA1loop1,:);
%         [AnnC1,AnnV1] = AnnOp(4,TestVect);
%         [AnnC2,AnnV2] = AnnOp(4,AnnV1);
%         [CreC1,CreV1] = CreaOp(2,AnnV2);
%         [CreC2,CreV2] = CreaOp(2,CreV1);
%        for FqA1loop3 = 1:length(Basis)
%         if CreV2(1,:) == Basis(FqA1loop3,:)
%              FqAP2 = C2^2 * C4^2 * Wavefunc(FqA1loop1,1) * Wavefunc(FqA1loop3,1) * AnnC1 * AnnC2 * CreC1 * CreC2;
%              FqAMat2(FqA1loop1,FqA1loop3) = FqAMat2(FqA1loop1,FqA1loop3) +FqAP2;
%         end
%        end
% 
%     %%%% For the rest of the terms
% 
%     TestVectFAp1 = Basis(FqA1loop1,:);TestVectFAp8 = Basis(FqA1loop1,:);
%     TestVectFAp2 = Basis(FqA1loop1,:);TestVectFAp9 = Basis(FqA1loop1,:);
%     TestVectFAp3 = Basis(FqA1loop1,:);TestVectFAp10 = Basis(FqA1loop1,:);
%     TestVectFAp4 = Basis(FqA1loop1,:);TestVectFAp11 = Basis(FqA1loop1,:);
%     TestVectFAp5 = Basis(FqA1loop1,:);TestVectFAp12 = Basis(FqA1loop1,:);
%     TestVectFAp6 = Basis(FqA1loop1,:);TestVectFAp13 = Basis(FqA1loop1,:);
%     TestVectFAp7 = Basis(FqA1loop1,:);TestVectFAp14 = Basis(FqA1loop1,:);
%                                       TestVectFAp15 = Basis(FqA1loop1,:);
% 
%     [P1Vect,P1Const] = CreaAnnDouble(TestVectFAp1,2,2,2,4);
%     for FqA1loop4 = 1:length(Basis)
%         if P1Vect == Basis(FqA1loop4,:)
%             FqMatAP1(FqA1loop1,FqA1loop4) = FqMatAP1(FqA1loop1,FqA1loop4) + ...
%                 (C2^3 * C4 * Wavefunc(FqA1loop1,1) * Wavefunc(FqA1loop4,1)*P1Const);
%         end
%     end
%        [P2Vect,P2Const] = CreaAnnDouble(TestVectFAp2,2,4,2,2);
%     for FqA1loop4 = 1:length(Basis)
%         if P2Vect == Basis(FqA1loop4,:)
%             FqMatAP1(FqA1loop1,FqA1loop4) = FqMatAP1(FqA1loop1,FqA1loop4) + ...
%                 (C2^3 * C4 * Wavefunc(FqA1loop1,1) * Wavefunc(FqA1loop4,1)*P2Const);
%         end
%     end
%     [P3Vect,P3Const] = CreaAnnDouble(TestVectFAp3,2,4,2,2);
%     for FqA1loop4 = 1:length(Basis)
%         if P3Vect == Basis(FqA1loop4,:)
%             FqMatAP1(FqA1loop1,FqA1loop4) = FqMatAP1(FqA1loop1,FqA1loop4) + ...
%                 (C2^3 * C4 * Wavefunc(FqA1loop1,1) * Wavefunc(FqA1loop4,1)*P3Const);
%         end
%     end
%     [P4Vect,P4Const] = CreaAnnDouble(TestVectFAp4,4,2,2,2);
%     for FqA1loop4 = 1:length(Basis)
%         if P4Vect == Basis(FqA1loop4,:)
%             FqMatAP1(FqA1loop1,FqA1loop4) = FqMatAP1(FqA1loop1,FqA1loop4) + ...
%                 (C2^3 * C4 * Wavefunc(FqA1loop1,1) * Wavefunc(FqA1loop4,1)*P4Const);
%         end
%     end
%     [P5Vect,P5Const] = CreaAnnDouble(TestVectFAp5,4,2,2,4);
%     for FqA1loop4 = 1:length(Basis)
%         if P5Vect == Basis(FqA1loop4,:)
%             FqMatAP1(FqA1loop1,FqA1loop4) = FqMatAP1(FqA1loop1,FqA1loop4) + ...
%                 (C2^3 * C4 * Wavefunc(FqA1loop1,1) * Wavefunc(FqA1loop4,1)*P5Const);
%         end
%     end
%     [P6Vect,P6Const] = CreaAnnDouble(TestVectFAp6,4,4,2,2);
%     for FqA1loop4 = 1:length(Basis)
%         if P6Vect == Basis(FqA1loop4,:)
%             FqMatAP1(FqA1loop1,FqA1loop4) = FqMatAP1(FqA1loop1,FqA1loop4) + ...
%                 (C2^3 * C4 * Wavefunc(FqA1loop1,1) * Wavefunc(FqA1loop4,1)*P6Const);
%         end
%     end
%     [P7Vect,P7Const] = CreaAnnDouble(TestVectFAp7,4,4,4,2);
%     for FqA1loop4 = 1:length(Basis)
%         if P7Vect == Basis(FqA1loop4,:)
%             FqMatAP1(FqA1loop1,FqA1loop4) = FqMatAP1(FqA1loop1,FqA1loop4) + ...
%                 (C2^3 * C4 * Wavefunc(FqA1loop1,1) * Wavefunc(FqA1loop4,1)*P7Const);
%         end
%     end
% 
%     [P8Vect,P8Const] = CreaAnnDouble(TestVectFAp8,2,4,4,4);
%     for FqA1loop4 = 1:length(Basis)
%         if P8Vect == Basis(FqA1loop4,:)
%             FqMatAP1(FqA1loop1,FqA1loop4) = FqMatAP1(FqA1loop1,FqA1loop4) + ...
%                 (C2^1 * C4^3 * Wavefunc(FqA1loop1,1) * Wavefunc(FqA1loop4,1)*P8Const);
%         end
%     end
%     [P9Vect,P9Const] = CreaAnnDouble(TestVectFAp9,4,2,2,4);
%     for FqA1loop4 = 1:length(Basis)
%         if P9Vect == Basis(FqA1loop4,:)
%             FqMatAP1(FqA1loop1,FqA1loop4) = FqMatAP1(FqA1loop1,FqA1loop4) + ...
%                 (C2^1 * C4^3 * Wavefunc(FqA1loop1,1) * Wavefunc(FqA1loop4,1)*P9Const);
%         end
%     end
%     [P10Vect,P10Const] = CreaAnnDouble(TestVectFAp10,4,4,4,2);
%     for FqA1loop4 = 1:length(Basis)
%         if P10Vect == Basis(FqA1loop4,:)
%             FqMatAP1(FqA1loop1,FqA1loop4) = FqMatAP1(FqA1loop1,FqA1loop4) + ...
%                 (C2^1 * C4^3 * Wavefunc(FqA1loop1,1) * Wavefunc(FqA1loop4,1)*P10Const);
%         end
%     end
% 
%     [P11Vect,P11Const] = CreaAnnSi(TestVectFAp11,2,4);
%     for FqA1loop4 = 1:length(Basis)
%         if P11Vect == Basis(FqA1loop4,:)
%             FqMatAP1(FqA1loop1,FqA1loop4) = FqMatAP1(FqA1loop1,FqA1loop4) + ...
%                 (C2^3 * C4^1 * Wavefunc(FqA1loop1,1) * Wavefunc(FqA1loop4,1)*P11Const);
%         end
%     end
%     [P12Vect,P12Const] = CreaAnnSi(TestVectFAp12,4,2);
%     for FqA1loop4 = 1:length(Basis)
%         if P12Vect == Basis(FqA1loop4,:)
%             FqMatAP1(FqA1loop1,FqA1loop4) = FqMatAP1(FqA1loop1,FqA1loop4) + ...
%                 (C2^3 * C4^1 * Wavefunc(FqA1loop1,1) * Wavefunc(FqA1loop4,1)*P12Const);
%         end
%     end
%     [P13Vect,P13Const] = CreaAnnSi(TestVectFAp13,4,4);
%     for FqA1loop4 = 1:length(Basis)
%         if P13Vect == Basis(FqA1loop4,:)
%             FqMatAP1(FqA1loop1,FqA1loop4) = FqMatAP1(FqA1loop1,FqA1loop4) + ...
%                 (C2^3 * C4^1 * Wavefunc(FqA1loop1,1) * Wavefunc(FqA1loop4,1)*P13Const);
%         end
%     end
%     [P14Vect,P14Const] = CreaAnnSi(TestVectFAp14,2,4);
%     for FqA1loop4 = 1:length(Basis)
%         if P14Vect == Basis(FqA1loop4,:)
%             FqMatAP1(FqA1loop1,FqA1loop4) = FqMatAP1(FqA1loop1,FqA1loop4) + ...
%                 (C2^1 * C4^3 * Wavefunc(FqA1loop1,1) * Wavefunc(FqA1loop4,1)*P14Const);
%         end
%     end
%     [P15Vect,P15Const] = CreaAnnSi(TestVectFAp15,4,2);
%     for FqA1loop4 = 1:length(Basis)
%         if P15Vect == Basis(FqA1loop4,:)
%             FqMatAP1(FqA1loop1,FqA1loop4) = FqMatAP1(FqA1loop1,FqA1loop4) + ...
%                 (C2^3 * C4^1 * Wavefunc(FqA1loop1,1) * Wavefunc(FqA1loop4,1)*P15Const);
%         end
%     end
% 
% end
% FqTotMat = FqAMat1+FqAMat2+FqMatAP1;
% FqA = sum(sum(FqAMat1+FqAMat2+FqMatAP1))
% %
% %% FqB
% %%%FqB1
% FqB1Mat = zeros(192);
% for FqB1loop1 = 1:length(Basis)
%     FqB1ann = zeros(1,22); FqB1ann(1,a2) = 1;
%     FqB1cre = zeros(1,22); FqB1cre(1,a2) = 1;
% 
%     TestVect = Basis(FqB1loop1,:);
%     annTestVect = TestVect - FqB1ann;
%     if any(annTestVect == -1)
%         continue
%     else creaannVect = annTestVect + FqB1cre;
%     end
%     for FqB1loop2 = 1:length(Basis)
%         if creaannVect == Basis(FqB1loop2,:)
%             Val = sqrt(Basis(FqB1loop1,a2)) * sqrt(annTestVect(1,a2)+1) * C2^2 * Wavefunc(FqB1loop1,1) * Wavefunc(FqB1loop2,1);
%             FqB1Mat(FqB1loop1,FqB1loop2) = FqB1Mat(FqB1loop1,FqB1loop2) + Val;
%         end
%     end
% end
% %%%FqB2
% FqB2Mat = zeros(192);
% for FqB2loop1 = 1:length(Basis)
%     FqB2ann = zeros(1,22); FqB2ann(1,a4) = 1;
%     FqB2cre = zeros(1,22); FqB2cre(1,a2) = 1;
% 
%     TestVect = Basis(FqB2loop1,:);
%     annTestVect = TestVect - FqB2ann;
%     if any(annTestVect == -1)
%         continue
%     else creaannVect = annTestVect + FqB2cre;
%     end
%     for FqB2loop2 = 1:length(Basis)
%         if creaannVect == Basis(FqB2loop2,:)
%             Val = sqrt(Basis(FqB2loop1,a4)) * sqrt(annTestVect(1,a2)+1) * C2*C4 * Wavefunc(FqB1loop1,1) * Wavefunc(FqB1loop2,1);
%             FqB2Mat(FqB2loop1,FqB2loop2) = FqB2Mat(FqB2loop1,FqB2loop2) + Val;
%         end
%     end
% end
% %%%FqB3
% FqB3Mat = zeros(192);
% for FqB3loop1 = 1:length(Basis)
%     FqB3ann = zeros(1,22); FqB3ann(1,a2) = 1;
%     FqB3cre = zeros(1,22); FqB3cre(1,a4) = 1;
% 
%     TestVect = Basis(FqB3loop1,:);
%     annTestVect = TestVect - FqB3ann;
%     if any(annTestVect == -1)
%         continue
%     else creaannVect = annTestVect + FqB3cre;
%     end
%     for FqB3loop2 = 1:length(Basis)
%         if creaannVect == Basis(FqB3loop2,:)
%             Val = sqrt(Basis(FqB3loop1,a2)) * sqrt(annTestVect(1,a4)+1) * C2*C4 * Wavefunc(FqB1loop1,1) * Wavefunc(FqB1loop2,1);
%             FqB3Mat(FqB3loop1,FqB3loop2) = FqB1Mat(FqB3loop1,FqB3loop2) + Val;
%         end
%     end
% end
% %%%FqB4
% FqB4Mat = zeros(192);
% for FqB4loop1 = 1:length(Basis)
%     FqB4ann = zeros(1,22); FqB4ann(1,a4) = 1;
%     FqB4cre = zeros(1,22); FqB4cre(1,a4) = 1;
% 
%     TestVect = Basis(FqB4loop1,:);
%     annTestVect = TestVect - FqB4ann;
%     if any(annTestVect == -1)
%         continue
%     else creaannVect = annTestVect + FqB4cre;
%     end
%     for FqB4loop2 = 1:length(Basis)
%         if creaannVect == Basis(FqB4loop2,:)
%             Val = sqrt(Basis(FqB4loop1,a4)) * sqrt(annTestVect(1,a4)+1) * C4^2 * Wavefunc(FqB1loop1,1) * Wavefunc(FqB1loop2,1);
%             FqB4Mat(FqB4loop1,FqB4loop2) = FqB4Mat(FqB4loop1,FqB4loop2) + Val;
%         end
%     end
% end
% FqBMatTot = FqB1Mat+FqB2Mat+FqB3Mat+FqB4Mat;
% FqB = sum(sum(FqBMatTot.^2));
% 
% %
% %
% %
% %
% %
% %
% % %%
% %
% Fq = 4*(FqA-FqB)
% TempMat(1,OmegaC) = Omega;
% TempMat(2,OmegaC) = Fq;










toc