%%%%%% Section selects the number of states to calc depending on if they
%%%%%% exceed a pre-defined tolerance, ie > 0.97



clc; clear all;
tic

KPos = [ 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22; %Matlab position
    0 0 0 0 0 0 0 0 0 0  0  1  1  1  1  1  1  1  1  1  1  1 ; %n value
    -1 0 1 2 3 4 5 6 7 8  9 -1  0  1  2  3  4  5  6  7  8  9]; %Mt value
format long;
A = 0.03 ;
g = 0.4;

Nt = 6;
%%%% Omega = 0.9375722836016 %best so far for g = 0.4
%%%% Omega = 0.82856708053814 for g = 1
OmegaC = 0;
for Omega = 0.9375722836016
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
GsVect = eigVectSPDM(:,maxEigPos)
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


Ex1VectTemp(:,1) = Ex1Vect(:,1);
Ex1VectTemp(:,2) = abs(Ex1Vect(:,1));
[Ex1Val,Ex1Pos] = max(Ex1VectTemp(:,2));
Ex1Val = Ex1VectTemp(Ex1Pos,1);
Ex1VectTemp(Ex1Pos,2) = 0;
[Ex2Val,Ex2Pos]  =max(Ex1VectTemp(:,2));
Ex2Val = Ex1VectTemp(Ex2Pos,1);
Ex1VectTemp(Ex2Pos,2) = 0;
[Ex3Val,Ex3Pos] = max(Ex1VectTemp(:,2));
Ex3Val = Ex1VectTemp(Ex3Pos,1);



%%%%% USEFUL FOR TROUBLESHOOTING BUT A BIT OF A RELIC IF THE SELECTION TOO
%%%%% LWOEKS
% Vects = [GsVect Ex1Vect]
% disp('Gs Val and Pos')
% GsData = [Gs1Val Gs1Pos;Gs2Val Gs2Pos;Gs3Val Gs3Pos]
% GsFid = (Gs1Val.^2)+(Gs2Val.^2)+(Gs3Val.^2)
% disp('Ex1 Val and Pos')
% Ex1Data = [Ex1Val  Ex1Pos; Ex2Val  Ex2Pos; Ex3Val  Ex3Pos]
% Ex1Fid = Ex1Val.^2 + Ex2Val.^2 + Ex3Val.^2

Tol = 0.995; %Tolerence to beat
%%%%GS selection tool
if Gs1Val.^2 > Tol
    Gs1  = 1; Gs2 = 0; Gs3 = 0;
elseif Gs1Val.^2 + Gs2Val.^2 > Tol
    Gs1 = 1; Gs2 = 1; Gs3 = 0;
else Gs1 =1;Gs2 =1;Gs3 = 1;
end

if Ex1Val.^2 > Tol
    Ex1  = 1; Ex2 = 0; Ex3 = 0;
elseif Ex1Val.^2 + Ex2Val.^2 > Tol
    Ex1 = 1; Ex2 = 1; Ex3 = 0;
else Ex1 =1;Ex2 =1;Ex3 = 1;
end
%Gs1
%Gs2
%Gs3
%Ex1
%Ex2
%Ex3


%% TrixTri FUNCTION
if (Gs1 == 1) && (Gs2 == 1) && (Gs3 == 1) && (Ex1 == 1) && (Ex2 == 1) && (Ex3 == 1)
    [Created60Kets,Created42Kets,Created24Kets,Created06Kets] =...
        FnTrixTri(Gs1Val,Gs2Val,Gs3Val,Gs1Pos,Gs2Pos,Gs3Pos,Ex1Val,Ex2Val,Ex3Val,Ex1Pos,Ex2Pos,Ex3Pos);
%    disp('TrixTri');
end
%% TrixBi FUNCTION
if (Gs1 == 1) && (Gs2 == 1) && (Gs3 == 1) && (Ex1 == 1) && (Ex2 == 1) && (Ex3 == 0)
    [Created60Kets,Created42Kets,Created24Kets,Created06Kets] =...
        FnTrixBi(Gs1Val,Gs2Val,Gs3Val,Gs1Pos,Gs2Pos,Gs3Pos,Ex1Val,Ex2Val,Ex1Pos,Ex2Pos);
%    disp('TrixBi');
end
%% TrixSi
if (Gs1 == 1) && (Gs2 == 1) && (Gs3 == 1) && (Ex1 == 1) && (Ex2 == 0) && (Ex3 == 0)
    [Created60Kets,Created42Kets,Created24Kets,Created06Kets] =...
        FnTrixSi(Gs1Val,Gs2Val,Gs3Val,Gs1Pos,Gs2Pos,Gs3Pos,Ex1Val,Ex1Pos);
%    disp('TrixSi');
end
%% BixTri
if (Gs1 == 1) && (Gs2 == 1) && (Gs3 == 0) && (Ex1 == 1) && (Ex2 == 1) && (Ex3 == 1)
    [Created60Kets,Created42Kets,Created24Kets,Created06Kets] =...
        FnBixTri(Gs1Val,Gs2Val,Gs1Pos,Gs2Pos,Ex1Val,Ex2Val,Ex3Val,Ex1Pos,Ex2Pos,Ex3Pos);
%    disp('BixTri');
end
%% BixBi
if (Gs1 == 1) && (Gs2 == 1) && (Gs3 == 0) && (Ex1 == 1) && (Ex2 == 1) && (Ex3 == 0)
    [Created60Kets,Created42Kets,Created24Kets,Created06Kets] =...
        FnBixBi(Gs1Val,Gs2Val,Gs1Pos,Gs2Pos,Ex1Val,Ex2Val,Ex1Pos,Ex2Pos);
%    disp('BixBi');
end
%% BixSi
if (Gs1 == 1) && (Gs2 == 1) && (Gs3 == 0) && (Ex1 == 1) && (Ex2 == 0) && (Ex3 == 0)
    disp('MISSING BixSi')
    
end
%% SixBi
if (Gs1 == 1) && (Gs2 == 0) && (Gs3 == 0) && (Ex1 == 1) && (Ex2 == 1) && (Ex3 == 0)
    disp('MISSING SixBi')
    
end
%% SixSi
if (Gs1 == 1) && (Gs2 == 0) && (Gs3 == 0) && (Ex1 == 1) && (Ex2 == 0) && (Ex3 == 0)
    disp('MISSING SixSi')
    
end
%% SixTri
if (Gs1 == 1) && (Gs2 == 0) && (Gs3 == 0) && (Ex1 == 1) && (Ex2 == 1) && (Ex3 == 1)
    disp('MISSING BixT')
    
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


format long
prob60 = (sum(mat60(:,1)).^2);
prob42 = (sum(mat42(:,1)).^2);
prob24 = (sum(mat24(:,1)).^2);
prob06 = (sum(mat06(:,1)).^2);
Fid = prob60 + prob42 + prob24 + prob06

bar([prob60  prob42  prob24  prob06]);
%txt = ['g = ' num2str(g)];txt2 = ['Fidelity = ' num2str(Fid)];
%text(2,0.3,txt)
%text(2,0.25,txt2)
dim = [0.4 .5 .3 .3];
str = ['g = ' num2str(g) char(10) 'A= ' num2str(A) char(10)... 
    'Omega = ' num2str(Omega) char(10) 'Fidelity = ' num2str(Fid)]; % char(10) starts a new line
annotation('textbox',dim,'String',str, 'EdgeColor','none')

Mat(1,OmegaC) = Omega;
Mat(2,OmegaC) = prob60;
Mat(3,OmegaC) = prob42;
Mat(4,OmegaC) = prob24;
Mat(5,OmegaC) = prob06;
Mat(6,OmegaC) = Fid;

end












toc