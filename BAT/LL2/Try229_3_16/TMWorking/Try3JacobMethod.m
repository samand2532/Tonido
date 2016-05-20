%%%%%% Tri, Bi expansion, new method for Tri-Bi expansion


clc; clear all;
tic

KPos = [ 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22; %Matlab position
    0 0 0 0 0 0 0 0 0 0  0  1  1  1  1  1  1  1  1  1  1  1 ; %n value
    -1 0 1 2 3 4 5 6 7 8  9 -1  0  1  2  3  4  5  6  7  8  9]; %Mt value

A = 0.03;
g = 0.4;

Nt = 6;

OmegaLoop = 0;
%for Omega = 0.93757:0.00000001:0.937575
%Omega = 0.9389% GOOD FOR g = 0.4
Omega = 0.94
OmegaLoop = OmegaLoop + 1;



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
format short

[eigVectSPDM,eigValSPDM] = eig(SPDMmat);

%%%%% finds 2 largest (+ve) eigenstates
%%% by diaganolising (MATLAB DEFINITION) the eigValSPDM mat, finding the
%%% corresponding row, then overwriting it with 0 then
%%% finding the next max.

diaEigVal = diag(eigValSPDM);
[maxEigVal,maxEigPos] = max(diaEigVal);
diaEigVal(maxEigPos,1) = 0;
[SecmaxEigVal,SecmaxEigPos] = max(diaEigVal);
GsVect = eigVectSPDM(:,maxEigPos);
Ex1Vect = eigVectSPDM(:,SecmaxEigPos);
GsVectTemp(:,1) = (GsVect(:));
GsVectTemp(:,2) = abs(GsVect(:));
[Gs1Val,Gs1Pos] = max(GsVectTemp(:,2));
Gs1Val = GsVectTemp(Gs1Pos,1);
GsVectTemp(Gs1Pos,:) = 0;
[Gs2Val,Gs2Pos] = max(GsVectTemp(:,2));
Gs2Val = GsVectTemp(Gs2Pos,1);
GsVectTemp(Gs2Pos,:) = 0;
[Gs3Val,Gs3Pos] = max(GsVectTemp(:,2));
Gs3Val = GsVectTemp(Gs3Pos,1);
GsVectTemp(Gs3Pos,:) = 0;
Vects = [GsVect Ex1Vect]
disp('Gs Val and Pos')
GsData = [Gs1Val Gs1Pos;Gs2Val Gs2Pos;Gs3Val Gs3Pos]
GsFid = (Gs1Val.^2)+(Gs2Val.^2)+(Gs3Val.^2)

Ex1VectTemp(:,1) = Ex1Vect(:,1);
Ex1VectTemp(:,2) = abs(Ex1Vect(:,1));
[Ex1Val,Ex1Pos] = max(Ex1VectTemp(:,2));
Ex1Val = Ex1VectTemp(Ex1Pos,1);
Ex1VectTemp(Ex1Pos,2) = 0;
[Ex2Val,Ex2Pos] = max(Ex1VectTemp(:,2));
disp('Ex1 Val and Pos')
Ex1Data = [Ex1Val  Ex1Pos; Ex2Val  Ex2Pos]
Ex1Fid = Ex1Val.^2 + Ex2Val.^2

TEMPMAT(OmegaLoop,1) = Omega;
TEMPMAT(OmegaLoop,2) = Gs1Pos;
TEMPMAT(OmegaLoop,3) = Gs2Pos;
TEMPMAT(OmegaLoop,4) = Gs3Pos;

AllCreatedKets = zeros(510,23);
Created60Kets = zeros(210,23);
Created42Kets = zeros(210,23);
Created24Kets = zeros(210,23);
Created06Kets = zeros(210,23);
aa=0;bb=0;cc=0;dd=0;ee=0;ff=0;
for loopi = 0:6
    for loopj = 0:6
        for loopk = 0:6
            for loopp = 0:6
                for loopq = 0:6
                    
                    if (loopi+loopj+loopk+loopp+loopq) == 6
                        aa=aa+1;
                        N1 = (loopi+loopj+loopk);
                        N2 = Nt - N1;
                        fN1 = factorial(N1);
                        fN2 = factorial(N2);
                        fi = factorial(loopi);
                        fj = factorial(loopj);
                        fk = factorial(loopk);
                        fp = factorial(loopp);
                        fq = factorial(loopq);
                        
                        TriBi(aa,1) = N1;
                        TriBi(aa,2) = N2;
                        TriBi(aa,3) = loopi;
                        TriBi(aa,4) = loopj;
                        TriBi(aa,5) = loopk;
                        TriBi(aa,6) = loopp;
                        TriBi(aa,7) = loopq;
                        
                        AllCreatedKets(aa,Gs1Pos+1) = loopi;
                        AllCreatedKets(aa,Gs2Pos+1) = loopj;
                        AllCreatedKets(aa,Gs3Pos+1) = loopk;
                        AllCreatedKets(aa,Ex1Pos+1) = loopp;
                        AllCreatedKets(aa,Ex2Pos+1) = loopq;
                        
                        AllCreatedKets(aa,1) = (1/(sqrt(fN1*fN2))) * (fN1 / (fi*fj*fk))...
                            * (Gs1Val^loopi)*(Gs2Val^loopj)*(Gs3Val^loopk)*sqrt(fi)*sqrt(fj)*sqrt(fk)...
                            * nchoosek(N2,loopp)*(Ex1Val^loopp)*(Ex2Val^loopq)*sqrt(fp)*sqrt(fq);
                           
                        if N1 == 6 && N2 ==0
                            bb=bb+1;
                            Created60Kets(bb,Gs1Pos+1) = loopi;
                            Created60Kets(bb,Gs2Pos+1) = loopj;
                            Created60Kets(bb,Gs3Pos+1) = loopk;
                            Created60Kets(bb,Ex1Pos+1) = loopp;
                            Created60Kets(bb,Ex2Pos+1) = loopq;
                            
                            Created60Kets(bb,1) = (1/(sqrt(fN1*fN2))) * (fN1 / (fi*fj*fk))...
                                * (Gs1Val^loopi)*(Gs2Val^loopj)*(Gs3Val^loopk)*sqrt(fi)*sqrt(fj)*sqrt(fk)...
                                * nchoosek(N2,loopp)*(Ex1Val^loopp)*(Ex2Val^loopq)*sqrt(fp)*sqrt(fq);
                        elseif N1 == 4 && N2 ==2
                            cc=cc+1;
                            Created42Kets(cc,Gs1Pos+1) = loopi;
                            Created42Kets(cc,Gs2Pos+1) = loopj;
                            Created42Kets(cc,Gs3Pos+1) = loopk;
                            Created42Kets(cc,Ex1Pos+1) = loopp;
                            Created42Kets(cc,Ex2Pos+1) = loopq;
                            
                            Created42Kets(cc,1) = (1/(sqrt(fN1*fN2))) * (fN1 / (fi*fj*fk))...
                                * (Gs1Val^loopi)*(Gs2Val^loopj)*(Gs3Val^loopk)*sqrt(fi)*sqrt(fj)*sqrt(fk)...
                                * nchoosek(N2,loopp)*(Ex1Val^loopp)*(Ex2Val^loopq)*sqrt(fp)*sqrt(fq);
                        elseif N1 == 2 && N2 ==4
                            dd=dd+1;
                            Created24Kets(dd,Gs1Pos+1) = loopi;
                            Created24Kets(dd,Gs2Pos+1) = loopj;
                            Created24Kets(dd,Gs3Pos+1) = loopk;
                            Created24Kets(dd,Ex1Pos+1) = loopp;
                            Created24Kets(dd,Ex2Pos+1) = loopq;
                            
                            Created24Kets(dd,1) = (1/(sqrt(fN1*fN2))) * (fN1 / (fi*fj*fk))...
                                * (Gs1Val^loopi)*(Gs2Val^loopj)*(Gs3Val^loopk)*sqrt(fi)*sqrt(fj)*sqrt(fk)...
                                * nchoosek(N2,loopp)*(Ex1Val^loopp)*(Ex2Val^loopq)*sqrt(fp)*sqrt(fq);
                        elseif N1 == 0 && N2 ==6
                            ee=ee+1;
                            Created06Kets(ee,Gs1Pos+1) = loopi;
                            Created06Kets(ee,Gs2Pos+1) = loopj;
                            Created06Kets(ee,Gs3Pos+1) = loopk;
                            Created06Kets(ee,Ex1Pos+1) = loopp;
                            Created06Kets(ee,Ex2Pos+1) = loopq;
                            
                            Created06Kets(ee,1) = (1/(sqrt(fN1*fN2))) * (fN1 / (fi*fj*fk))...
                                * (Gs1Val^loopi)*(Gs2Val^loopj)*(Gs3Val^loopk)*sqrt(fi)*sqrt(fj)*sqrt(fk)...
                                * nchoosek(N2,loopq)*(Ex1Val^loopq)*(Ex2Val^loopp)*sqrt(fp)*sqrt(fq);   
                        end
                    end
                end
            end
        end
    end
end
wave = [OrigCoeff Basis];
[x60,y60] = size(Created60Kets);
[x42,y42] = size(Created42Kets);
[x24,y24] = size(Created24Kets);
[x06,y06] = size(Created06Kets);
xx=0;
for CountA = 1:x60
   for Basis = 1:length(wave)
    if wave(Basis,2:end) == Created60Kets(CountA,2:end)
        xx=xx+1;
        mat60(xx,1) = wave(Basis,1)*Created60Kets(CountA,1);
    end
   end
end
xx=0;
for CountA = 1:x42
   for Basis = 1:length(wave)
    if wave(Basis,2:end) == Created42Kets(CountA,2:end)
        xx=xx+1;
        mat42(xx,1) = wave(Basis,1)*Created42Kets(CountA,1);
    end
   end
end
xx=0;
for CountA = 1:x24
   for Basis = 1:length(wave)
    if wave(Basis,2:end) == Created24Kets(CountA,2:end)
        xx=xx+1;
        mat24(xx,1) = wave(Basis,1)*Created24Kets(CountA,1);
    end
   end
end
xx=0;
for CountA = 1:x06
   for Basis = 1:length(wave)
    if wave(Basis,2:end) == Created06Kets(CountA,2:end)
        xx=xx+1;
        mat06(xx,1) = wave(Basis,1)*Created06Kets(CountA,1);
    end
   end
end



prob60 = sum(mat60(:,1).^2);
prob42 = sum(mat42(:,1).^2);
prob24 = sum(mat24(:,1).^2);
prob06 = sum(mat06(:,1).^2);



bar([prob60  prob42  prob24  prob06])
Fid = prob60 + prob42 + prob24 + prob06










toc