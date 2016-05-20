clc; clear all;
OmCount = 0;
OmMat = zeros(16300,8);
tic

    KPos = [ 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19; %Matlab position
        0 0 1 0 1 0 1 0 1 0  1  0  1  0  1  0  1  0  1 ; %n value
        -1 0 0 1 1 2 2 3 3 4  4  5  5  6  6  7  7  8  8]; %Mt value
    
    MomValue = [ -1 0 0 1 1 2 2 3 3 4 4 5 5 6 6 7 7 8 8];
    
    g = 1; N = 6; A = 0.015; LLVal = 1; % 1 = LLL, 2 = LL2
    format long
    
    if LLVal == 1
     for OmegaVal = 0.776755;%0.6:0.01:0.85%0.7771142:0.00000001:0.7772192;
        OmCount = OmCount+1
        
        AllBras = csvread('StatesLLL.csv');
        AllKets = AllBras';
        
        
        NMat = eye(length(AllBras)).*N;
        for BraL = 1:length(AllBras)
            for KetL = 1:length(AllKets)
                InitialBraL = AllBras(BraL,:);
                InitialKetL = AllKets(:,KetL);
                
                if InitialBraL == InitialKetL'
                    LMat(BraL,KetL) = dot(InitialKetL, MomValue);
                end
                
            end
        end
        
        deltaINT = csvread('LagIntLLL.csv');
        deltaPOT = csvread('LagPotLLL.csv');
        UMatNogConst = zeros(length(AllKets));
        
        l1Vect = zeros(size(InitialBraL))';
        l2Vect = zeros(size(InitialBraL))';
        k1Vect = zeros(size(InitialBraL))';
        k2Vect = zeros(size(InitialBraL))';
        %%%% INT PART
        for BrasU = 1:length(AllBras)
            BrasU;
            for KetsU = 1:length(AllBras)
                InitialBraU = AllBras(BrasU,:);
                InitialKetU = AllKets(:,KetsU);
                
                for deltaCountINT = 1:length(deltaINT)
                    k1 = deltaINT(deltaCountINT,1);
                    k2 = deltaINT(deltaCountINT,2);
                    l1 = deltaINT(deltaCountINT,3);
                    l2 = deltaINT(deltaCountINT,4);
                    nk1 = deltaINT(deltaCountINT,5);
                    nk2 = deltaINT(deltaCountINT,6);
                    n11 = deltaINT(deltaCountINT,7);
                    nl2 = deltaINT(deltaCountINT,8);
                    mk1 = deltaINT(deltaCountINT,9);
                    mk2 = deltaINT(deltaCountINT,10);
                    ml1 = deltaINT(deltaCountINT,11);
                    ml2 = deltaINT(deltaCountINT,12);
                    
                    l1Vect = l1Vect.*0;
                    l2Vect = l2Vect.*0;
                    k1Vect = k1Vect.*0;
                    k2Vect = k2Vect.*0;
                    
                    l1Vect(l1,1) = -1;
                    l2Vect(l2,1) = -1;
                    k1Vect(k1,1) = 1;
                    k2Vect(k2,1) = 1;
                    
                    l2Const = sqrt(InitialKetU(l2,1));
                    l2trans  = InitialKetU + l2Vect;
                    l1Const = sqrt(l2trans(l1,1));
                    l1trans = l2trans + l1Vect;
                    k2Const = sqrt(l1trans(k2,1)+1);
                    k2trans = l1trans + k2Vect;
                    k1Const = sqrt(k2trans(k1,1)+1);
                    FinalKet = k2trans + k1Vect;
                    
                    if FinalKet == InitialBraU'
                        
                        OpConst = l2Const*l1Const*k2Const*k1Const;
                        OtherConst = deltaINT(deltaCountINT,13);
                        UMatNogConst(BrasU,KetsU) = UMatNogConst(BrasU,KetsU) + (OtherConst*OpConst);
                    end
                end
            end
        end
        %%%%% POT TERM
        
        k1Vect = zeros(size(InitialBraL))';
        k2Vect = zeros(size(InitialBraL))';
        
        VMatNoA = zeros(length(AllKets));
        
        for BrasV = 1:length(AllBras)
            for KetsV = 1:length(AllBras)
                InitialBraV = AllBras(BrasV,:);
                InitialKetV = AllKets(:,KetsV);
                
                for deltaCountPOT = 1:length(deltaPOT)
                    
                    k1 = deltaPOT(deltaCountPOT,1);
                    k2 = deltaPOT(deltaCountPOT,2);
                    nk1 = deltaPOT(deltaCountPOT,3);
                    nk2 = deltaPOT(deltaCountPOT,4);
                    mk1 = deltaPOT(deltaCountPOT,5);
                    mk2 = deltaPOT(deltaCountPOT,6);
                    
                    k1Vect = k1Vect.*0;
                    k2Vect = k2Vect.*0;
                    
                    k1Vect(k1,1) = 1;
                    k2Vect(k2,1) = -1;
                    
                    k2Const = sqrt(InitialKetV(k2,1));
                    k2trans = InitialKetV + k2Vect;
                    k1Const = sqrt(k2trans(k1,1)+1);
                    FinalKetV = k2trans + k1Vect;
                    
                    if FinalKetV == InitialBraV'
                        
                        OpConst = k2Const*k1Const;
                        OtherConst = deltaPOT(deltaCountPOT,7);
                        %I1 = factorial((abs(mk1)+abs(mk2)+2)/2);
                        %PiT = sqrt(1/(factorial(abs(mk1))*factorial(abs(mk2))));
                        %OtherConst = I1*PiT;
                        VMatNoA(BrasV,KetsV) = VMatNoA(BrasV,KetsV) + (OtherConst*OpConst);
                    end
                end
            end
        end
        
        
        
        n = 0;
        for Omega = OmegaVal;
            n=n+1;
            
            NAll = NMat;
            LAll = LMat*(1-Omega);
            UAll = UMatNogConst.*(g/(4*pi));
            UAlla=round(UAll,8);
            VMat = VMatNoA.*A;
            %VMat = round(VMata,8);
            Total = NAll +UAlla+(LAll)+VMat;
            Eig = eig(Total);
            [V D] = eig(Total);
            
            
            OmegaIncreaseMatrix(1,n) = Omega;
            OmegaEnergyMatrix(:,n) = Eig;
            %%% Following line selects only the lowest 8 eigenstates
            OmegaMatRestricted = OmegaEnergyMatrix(1:9,:);
            
            
        end
        VTest=V;     
        
    
   
    
    AllBras = csvread('StatesLLL.csv');
    AllKets = AllBras';
    V = VTest;
    eigenstates = V(:,1);
    format long
    eigMat = zeros(9);
    
    creaVect = zeros(19,1);
    annVect = zeros(19,1);
    
    for Kets = 1:39;
        for k = [2 4 6 8 10 12 14 16 18];
            for l = [2 4 6 8 10 12 14 16 18];
                C1 = eigenstates(Kets,1);
                InitialKet = AllKets(:,Kets);
                creaVect = creaVect.*0;
                annVect = annVect.*0;
                creaVect(k,1) = 1;
                annVect(l,1) = -1;
                
                annConst = sqrt(InitialKet(l,1));
                annTrans = InitialKet + annVect;
                creaConst = sqrt(annTrans(k,1)+1);
                FinalKet = annTrans + creaVect;
                a=0;
                for a = 1:39
                    
                    Bra = AllBras(a,:);
                    C2 = eigenstates(a,1);
                    if Bra' == FinalKet
                        eigMat(k/2,l/2) = eigMat(k/2,l/2) + (creaConst*annConst*C1*C2);
                    end
                end
            end
        end
    end
    
    eigMat;
    [V1 D] = eig(eigMat);
    dia = diag(eigMat);
    N = 1/sqrt(sum(dia).^2);
    Norm = dia.*N;
    SumNormEigVal = sum(Norm );
    
    OccupationMat(OmCount,1) = OmegaVal;
    OccupationMat(OmCount,2) = dia(1,1)/N;
    OccupationMat(OmCount,3) = dia(2,1)/N;
    OccupationMat(OmCount,4) = (dia(1,1)+dia(2,1))/N;
    
    D;
    
    V1;
    
    [maxNum, maxIndex] = max(D(:));
    [row, col] = ind2sub(size(D),maxIndex);
    
    ABC = sort(V1(:,col));
    xABC = ABC(1,1);
    x1ABC = ABC(9,1);
    
    Nt = 6;
    
    %%%% CHANGE THESE VALUES ONLY!!!!%%%%
    x = xABC;
    x1 = x1ABC;
    Eig = V(:,1);
    %%%%%%%%%%%%%%%%%%%%%%
    
    k=0;
    KetN0=zeros(1,10);
    for N0 = 0
        for p = 0:N0
            k=k+1;
            
            Nk = nchoosek(N0,p);
            Norm = 1/((sqrt(factorial(N0))).*(sqrt(factorial(Nt-N0))));
            
            
            KetN0(k,1) = Nk.*Norm.*(x^p).*(x1^(N0-p)).*(sqrt(factorial(p)))...
                .*sqrt(factorial(N0-p)).*sqrt(factorial(Nt-N0));
            KetN0(k,2) = p;
            KetN0(k,3) = 6-N0;
            KetN0(k,4) = N0-p;
        end
    end
    
    k=0;
    % KetN1 = zeros(2,10);
    % for N1 = 1
    %     for p = 0:N1
    %         k=k+1;
    %
    %          Nk = nchoosek(N1,p)
    %          SumTerm = Nk.*((x)^(N1-p)).*((x1).^p)
    %          Norm = 1/((sqrt(factorial(N1))).*(sqrt(factorial(6-N1))))
    %          TotConst = SumTerm.*Norm.*Nk
    %
    %
    %         KetN1(k,1) = Nk.*SumTerm.*Norm*(sqrt(factorial(p))*sqrt(factorial(6-N1))*sqrt(factorial(N1-p)));
    %         KetN1(k,2) = p;
    %         KetN1(k,3) = 6-N1;
    %         KetN1(k,4) = N1-p;
    %     end
    % end
    
    k=0;
    KetN2 = zeros(3,10);
    for N2 = 2
        for p = 0:N2
            k=k+1;
            
            Nk = nchoosek(N2,p);
            Norm = 1/((sqrt(factorial(N2))).*(sqrt(factorial(Nt-N2))));
            
            
            KetN2(k,1) = Nk.*Norm.*(x^p).*(x1^(N2-p)).*(sqrt(factorial(p)))...
                .*sqrt(factorial(N2-p)).*sqrt(factorial(Nt-N2));
            KetN2(k,2) = p;
            KetN2(k,3) = 6-N2;
            KetN2(k,4) = N2-p;
        end
    end
    
    
    k=0;
    % KetN3 = zeros(4,10);
    % for N3 = 3
    %     for p = 0:N3
    %         k=k+1;
    %
    %          Nk = nchoosek(N3,p)
    %          SumTerm = Nk.*((x)^(N3-p)).*((x1).^p)
    %          Norm = 1/((sqrt(factorial(N3))).*(sqrt(factorial(6-N3))))
    %          TotConst = SumTerm.*Norm.*Nk
    %
    %
    %         KetN3(k,1) = Nk.*SumTerm.*Norm*(sqrt(factorial(p))*sqrt(factorial(6-N3))*sqrt(factorial(N3-p)));
    %         KetN3(k,2) = p;
    %         KetN3(k,3) = 6-N3;
    %         KetN3(k,4) = N3-p;
    %     end
    % end
    
    k=0;
    KetN4 = zeros(5,10);
    for N4 = 4
        for p = 0:N4
            k=k+1;
            
            Nk = nchoosek(N4,p);
            Norm = 1/((sqrt(factorial(N4))).*(sqrt(factorial(Nt-N4))));
            
            
            KetN4(k,1) = Nk.*Norm.*(x^p).*(x1^(N4-p)).*(sqrt(factorial(p)))...
                .*sqrt(factorial(N4-p)).*sqrt(factorial(Nt-N4));
            KetN4(k,2) = p;
            KetN4(k,3) = 6-N4;
            KetN4(k,4) = N4-p;
        end
    end
    
    k=0;
    % KetN5 = zeros(6,10);
    % for N5 = 5
    %     for p = 0:N5
    %         k=k+1;
    %
    %          Nk = nchoosek(N5,p)
    %          SumTerm = Nk.*((x)^(N5-p)).*((x1).^p)
    %          Norm = 1/((sqrt(factorial(N5))).*(sqrt(factorial(6-N5))))
    %          TotConst = SumTerm.*Norm.*Nk
    %
    %
    %         KetN5(k,1) = Nk.*SumTerm.*Norm*(sqrt(factorial(p))*sqrt(factorial(6-N5))*sqrt(factorial(N5-p)));
    %         KetN5(k,2) = p;
    %         KetN5(k,3) = 6-N5;
    %         KetN5(k,4) = N5-p;
    %     end
    % end
    
    
    k=0;
    KetN6 = zeros(7,10);
    for N6 = 6
        for p = 0:N6
            k=k+1;
            
            Nk = nchoosek(N6,p);
            Norm = 1/((sqrt(factorial(N6))).*(sqrt(factorial(Nt-N6))));
            
            
            KetN6(k,1) = Nk.*Norm.*(x^p).*(x1^(N6-p)).*(sqrt(factorial(p)))...
                .*sqrt(factorial(N6-p)).*sqrt(factorial(Nt-N6));
            KetN6(k,2) = p;
            KetN6(k,3) = 6-N6;
            KetN6(k,4) = N6-p;
        end
    end
    
    BasisChangeKet = [KetN0;KetN2;KetN4;KetN6];
    
    %%%%%Above generates the basis transformation from |N_1,N_2> to the
    %%%%%momentum basis used initially, next part finds fidelity with riginal
    %%%%%wavefunction.
    
    
    OrigStates = csvread('StatesLLLnoN.csv');
    OrigWavefunction = [Eig(:,1) OrigStates];
    
    %%%%%% Ket0 Part
    l=0;m=0;
    KetN0SaveMat = zeros(10,1);
    for CountA = 1:length(OrigWavefunction)
        l=l+1;
        if KetN0(1,2:10) == OrigWavefunction(l,2:10)
            m=m+1;
            KetN0SaveMat(m,1) = KetN0(1,1)*OrigWavefunction(l,1) ;
        end
    end
    
    %%%%%%% Ket2 Part
    l=0;m=0;
    KetN2SaveMat = zeros(10,1);
    for CountA = 1:length(OrigWavefunction)
        l=l+1;
        if KetN2(1,2:10) == OrigWavefunction(l,2:10)
            m=m+1;
            KetN2SaveMat(m,1) = KetN2(1,1)*OrigWavefunction(l,1);
        elseif KetN2(2,2:10) == OrigWavefunction(l,2:10)
            m=m+1;
            KetN2SaveMat(m,1) = KetN2(2,1)*OrigWavefunction(l,1);
        elseif KetN2(3,2:10) == OrigWavefunction(l,2:10)
            m=m+1;
            KetN2SaveMat(m,1) = KetN2(3,1)*OrigWavefunction(l,1);
        end
    end
    
    %%%%%%% Ket4 Part
    l=0;m=0;
    KetN4SaveMat = zeros(10,1);
    for CountA = 1:length(OrigWavefunction)
        l=l+1;
        if KetN4(1,2:10) == OrigWavefunction(l,2:10)
            m=m+1;
            KetN4SaveMat(m,1) = KetN4(1,1)*OrigWavefunction(l,1);
        elseif KetN4(2,2:10) == OrigWavefunction(l,2:10)
            m=m+1;
            KetN4SaveMat(m,1) = KetN4(2,1)*OrigWavefunction(l,1);
        elseif KetN4(3,2:10) == OrigWavefunction(l,2:10)
            m=m+1;
            KetN4SaveMat(m,1) = KetN4(3,1)*OrigWavefunction(l,1);
        elseif KetN4(4,2:10) == OrigWavefunction(l,2:10)
            m=m+1;
            KetN4SaveMat(m,1) = KetN4(4,1)*OrigWavefunction(l,1);
        elseif KetN4(5,2:10) == OrigWavefunction(l,2:10)
            m=m+1;
            KetN4SaveMat(m,1) = KetN4(5,1)*OrigWavefunction(l,1);
        end
    end
    
    %%%%%Ket6 Part
    l=0;m=0;
    KetN6SaveMat = zeros(10,1);
    for CountA = 1:length(OrigWavefunction)
        l=l+1;
        if KetN6(1,2:10) == OrigWavefunction(l,2:10)
            m=m+1;
            KetN6(1,1);
            OrigWavefunction(l,1);
            KetN6SaveMat(m,1) = KetN6(1,1)*OrigWavefunction(l,1);
        elseif KetN6(2,2:10) == OrigWavefunction(l,2:10)
            m=m+1;
            KetN6(2,1);
            OrigWavefunction(l,1);
            KetN6SaveMat(m,1) = KetN6(2,1)*OrigWavefunction(l,1);
        elseif KetN6(3,2:10) == OrigWavefunction(l,2:10)
            m=m+1;
            KetN6(3,1);
            OrigWavefunction(l,1);
            KetN6SaveMat(m,1) = KetN6(3,1)*OrigWavefunction(l,1);
        elseif KetN6(4,2:10) == OrigWavefunction(l,2:10)
            m=m+1;
            KetN6(4,1);
            OrigWavefunction(l,1);
            KetN6SaveMat(m,1) = KetN6(4,1)*OrigWavefunction(l,1);
        elseif KetN6(5,2:10) == OrigWavefunction(l,2:10)
            m=m+1;
            KetN6(5,1);
            OrigWavefunction(l,1);
            KetN6SaveMat(m,1) = KetN6(5,1)*OrigWavefunction(l,1);
        elseif KetN6(6,2:10) == OrigWavefunction(l,2:10)
            m=m+1;
            KetN6(6,1);
            OrigWavefunction(l,1);
            KetN6SaveMat(m,1) = KetN6(6,1)*OrigWavefunction(l,1);
        elseif KetN6(7,2:10) == OrigWavefunction(l,2:10)
            m=m+1;
            KetN6(7,1);
            OrigWavefunction(l,1);
            KetN6SaveMat(m,1) = KetN6(7,1)*OrigWavefunction(l,1);
        end
    end
    
    KetN0Overlap = (abs(sum(KetN0SaveMat)))^2;
    KetN2Overlap = (abs(sum(KetN2SaveMat)))^2;
    KetN4Overlap = (abs(sum(KetN4SaveMat)))^2;
    KetN6Overlap = (abs(sum(KetN6SaveMat)))^2;
    
    Overlaps = [KetN0Overlap KetN2Overlap KetN4Overlap KetN6Overlap]
    bar(Overlaps)
    %%%%% OmVal|Overlap1|OL2|OL3|OL4|Fidelity|diff1;4|diff2;3|sumdiff
    OmMat(OmCount,1) = OmegaVal;
    OmMat(OmCount,2) = KetN0Overlap;
    OmMat(OmCount,3) = KetN2Overlap;
    OmMat(OmCount,4) = KetN4Overlap;
    OmMat(OmCount,5) = KetN6Overlap;
    OmMat(OmCount,6) = (KetN0Overlap+KetN2Overlap+KetN4Overlap+KetN6Overlap);
    OmMat(OmCount,7) = abs(KetN0Overlap - KetN6Overlap);
    OmMat(OmCount,8) = abs(KetN2Overlap - KetN4Overlap);
    OmMat(OmCount,9) = 0;
    OmMat(OmCount,10) = OmMat(OmCount,7) + OmMat(OmCount,8);
     end
    end
toc
