%%% Fresh Start of BatExpansion, using variable state expansions depending
%%% on if the state meets ( and exceeds a certian tolerence, ie 0.95)


clc; clear all;
tic
format long
KPos = [ 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22; %Matlab position
         0 0 0 0 0 0 0 0 0 0  0  1  1  1  1  1  1  1  1  1  1  1 ; %n value
        -1 0 1 2 3 4 5 6 7 8  9 -1  0  1  2  3  4  5  6  7  8  9]; %Mt value

A = 0.03;
g = 0.4;

Nt = 6;

OmegaLoop = 0;
%for Omega = 0.9:0.000001:0.96
%Omega = 0.9375% GOOD FOR g = 0.4
%Omega = 0.828 g=1 bat????

%Omega = 0.9371  %%%% GOOD for g=0.4 BAT
Omega = 0.9374779
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
Ex1VectTemp(Ex1Pos,:) = 0;
[Ex2Val,Ex2Pos] = max(Ex1VectTemp(:,2));
Ex2Val = Ex1VectTemp(Ex2Pos,1);
Ex1VectTemp(Ex2Pos,:) = 0;
[Ex3Val,Ex3Pos] = max(Ex1VectTemp(:,2));
disp('Ex1 Val and Pos')
Ex1Data = [Ex1Val  Ex1Pos; Ex2Val  Ex2Pos; Ex3Val  Ex3Pos]
Ex1Fid = Ex1Val.^2 + Ex2Val.^2 + Ex3Val.^2

% %%%% Trinomial-Trinomial expansion part
aa=0;
for gi = 0:6
    for gj = 0:6
        for gk = 0:6
            for ei = 0:6
                for ej = 0:6
                    for ek = 0:6
                        
                        if (gi+gj+gk+ei+ej+ek) == 6
                            aa=aa+1;
                            TriTriMat(aa,1) = (gi+gj+gk);
                            TriTriMat(aa,2) = (ei+ej+ek);
                            TriTriMat(aa,3) = gi;
                            TriTriMat(aa,4) = gj;
                            TriTriMat(aa,5) = gk;
                            TriTriMat(aa,6) = ei;
                            TriTriMat(aa,7) = ej;
                            TriTriMat(aa,8) = ek;
                            %TriTriMat of Form N1|N2|gi|gj|gk|ei|ej|ek
                        end
                    end
                end
            end
        end
    end
end

Created06Kets = zeros(100,23);
Created24Kets = zeros(100,23);
Created42Kets = zeros(100,23);
Created60Kets = zeros(100,23);
aa=0; bb=0;cc=0;dd=0;ee=0;
for CountA = 1:length(TriTriMat)
    N1 = TriTriMat(CountA,1);
    N2 = TriTriMat(CountA,2);
    gi = TriTriMat(CountA,3);
    gj = TriTriMat(CountA,4);
    gk = TriTriMat(CountA,5);
    ei = TriTriMat(CountA,6);
    ej = TriTriMat(CountA,7);
    ek = TriTriMat(CountA,8);

    fN1 = factorial(N1);
    fN2 = factorial(N2);
    fgi = factorial(gi);
    fgj = factorial(gj);
    fgk = factorial(gk);
    fei = factorial(ei);
    fej = factorial(ej);
    fek = factorial(ek);
    
    if N1 == 0 && N2 == 6
        aa=aa+1;    
        Created06Kets(aa,Gs1Pos+1) = gi;
        Created06Kets(aa,Gs2Pos+1) = gj;
        Created06Kets(aa,Gs3Pos+1) = gk;
        Created06Kets(aa,Ex1Pos+1) = ei;
        Created06Kets(aa,Ex2Pos+1) = ej;
        Created06Kets(aa,Ex3Pos+1) = ek;
        Created06Kets(aa,1) = (1/sqrt(fN1*fN2))*(Gs1Val^gi)*(Gs2Val^gj)*(Gs3Val^gk)...
            *(Ex1Val^ei)*(Ex2Val^ej)*(Ex3Val^ek)*(fN1/(fgi*fgj*fgk))*(fN2/(fei*fej*fek))...
            *sqrt(fgi)*sqrt(fgj)*sqrt(fgk)*sqrt(fei)*sqrt(fej)*sqrt(fek);
    end
    if N1 == 2 && N2 == 4
        bb=bb+1;    
        Created24Kets(bb,Gs1Pos+1) = gi;
        Created24Kets(bb,Gs2Pos+1) = gj;
        Created24Kets(bb,Gs3Pos+1) = gk;
        Created24Kets(bb,Ex1Pos+1) = ei;
        Created24Kets(bb,Ex2Pos+1) = ej;
        Created24Kets(bb,Ex3Pos+1) = ek;
        Created24Kets(bb,1) = (1/sqrt(fN1*fN2))*(Gs1Val^gi)*(Gs2Val^gj)*(Gs3Val^gk)...
            *(Ex1Val^ei)*(Ex2Val^ej)*(Ex3Val^ek)*(fN1/(fgi*fgj*fgk))*(fN2/(fei*fej*fek))...
            *sqrt(fgi)*sqrt(fgj)*sqrt(fgk)*sqrt(fei)*sqrt(fej)*sqrt(fek);
    end
    if N1 == 4 && N2 == 2
        cc=cc+1;    
        Created42Kets(cc,Gs1Pos+1) = gi;
        Created42Kets(cc,Gs2Pos+1) = gj;
        Created42Kets(cc,Gs3Pos+1) = gk;
        Created42Kets(cc,Ex1Pos+1) = ei;
        Created42Kets(cc,Ex2Pos+1) = ej;
        Created42Kets(cc,Ex3Pos+1) = ek;
        Created42Kets(cc,1) = (1/sqrt(fN1*fN2))*(Gs1Val^gi)*(Gs2Val^gj)*(Gs3Val^gk)...
            *(Ex1Val^ei)*(Ex2Val^ej)*(Ex3Val^ek)*(fN1/(fgi*fgj*fgk))*(fN2/(fei*fej*fek))...
            *sqrt(fgi)*sqrt(fgj)*sqrt(fgk)*sqrt(fei)*sqrt(fej)*sqrt(fek);
    end
    if N1 == 6 && N2 == 0
        dd=dd+1;    
        Created60Kets(dd,Gs1Pos+1) = gi;
        Created60Kets(dd,Gs2Pos+1) = gj;
        Created60Kets(dd,Gs3Pos+1) = gk;
        Created60Kets(dd,Ex1Pos+1) = ei;
        Created60Kets(dd,Ex2Pos+1) = ej;
        Created60Kets(dd,Ex3Pos+1) = ek;
        Created60Kets(dd,1) = (1/sqrt(fN1*fN2))*(Gs1Val^gi)*(Gs2Val^gj)*(Gs3Val^gk)...
            *(Ex1Val^ei)*(Ex2Val^ej)*(Ex3Val^ek)*(fN1/(fgi*fgj*fgk))*(fN2/(fei*fej*fek))...
            *sqrt(fgi)*sqrt(fgj)*sqrt(fgk)*sqrt(fei)*sqrt(fej)*sqrt(fek);
        
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
        WaveVal = wave(Basis,1)
        CreatedVal = Created60Kets(CountA,1)
        WavetimesCreated = wave(Basis,1)*Created60Kets(CountA,1)
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



prob60 = (sum(mat60(:,1)).^2)
prob42 = (sum(mat42(:,1)).^2)
prob24 = (sum(mat24(:,1)).^2)
prob06 = (sum(mat06(:,1)).^2)



bar([prob60  prob42  prob24  prob06])
Fid = prob60 + prob42 + prob24 + prob06












toc