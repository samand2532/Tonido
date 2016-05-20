clc; clear all;
APA = csvread('APA.csv');

AllBras = csvread('AllBrasLLL.csv');
AllKets = AllBras';
APA = csvread('APA.csv');
MomValue = [ -1 0 1 2 3 4 5 6 7 8 9 10 0]';
UMatNoConst = zeros(length(AllKets));

m1Vect = zeros(13,1);
m2Vect = zeros(13,1);
m3Vect = zeros(13,1);
m4Vect = zeros(13,1);

for BrasU = 1:length(AllBras)
    for KetsU = 1:length(AllKets)
        InitialBraU = AllBras(BrasU,:)';
        InitialKetU = AllKets(:,KetsU);
        
        for rr = 1:length(APA)
            m1 = APA(rr,1)+2;
            m2 = APA(rr,2)+2;
            m3 = APA(rr,3)+2;
            m4 = APA(rr,4)+2;
            m1Vect = m1Vect*0;
            m2Vect = m2Vect*0;
            m3Vect = m3Vect*0;
            m4Vect = m4Vect*0;
            m1Vect(m1,1) = 1;
            m2Vect(m2,1) = 1;
            m3Vect(m3,1) = -1;
            m4Vect(m4,1) = -1;
            
            m4Const = sqrt(InitialKetU(m4,1));
            m4trans = InitialKetU + m4Vect;
            m3Const = sqrt(m4trans(m3,1));
            m3trans = m4trans + m3Vect;
            m2Const = sqrt(m3trans(m2,1)+1);
            m2trans = m3trans + m2Vect;
            m1Const = sqrt(m2trans(m1,1)+1);
            FinalKetU = m2trans + m1Vect;
            
            if FinalKetU == InitialBraU
                OpConst = m1Const*m2Const*m3Const*m4Const
                I2 = gamma(((abs(m1-2))+(abs(m2-2))+(abs(m3-2))+(abs(m4-2)))/2 + 1)
                
                if APA(rr,5) == 0
                    nk1 = 0; nk2 = 0; nk3 = 0; nk4 = 0;
                elseif APA(rr,5) == 1
                    nk1 = 1; nk2 = 0; nk3 = 0; nk4 = 0;
                elseif APA(rr,5) == 2
                    nk1 = 0; nk2 = 1; nk3 = 0; nk4 = 0;
                elseif APA(rr,5) == 3
                    nk1 = 0; nk2 = 0; nk3 = 1; nk4 = 0;
                elseif APA(rr,5) == 4
                    nk1 = 0; nk2 = 0; nk3 = 0; nk4 = 1;
                end
                
                PiTerm = (factorial(nk1)*factorial(nk2)*factorial(nk3)*factorial(nk4))/((factorial(nk1 + abs(m1-2)))*(factorial(nk2 + abs(m2-2)))*(factorial(nk3 + abs(m3-2)))*(factorial(nk4 + abs(m4-2))))
                MomTerm = 1/(2^(dot(InitialKetU,abs(MomValue))/2))
                UMid = OpConst*I2*PiTerm*MomTerm;
                UMatNoConst(BrasU,KetsU) = UMatNoConst(BrasU,KetsU)+UMid;
                
            end
        end
    end
end

                UMatNoConst
                
                
                
                
                
                
                
                
                
                
                
                
                
                