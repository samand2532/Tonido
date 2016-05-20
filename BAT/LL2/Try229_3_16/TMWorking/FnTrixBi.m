function [ Created60Kets,Created42Kets,Created24Kets,Created06Kets ] =...
    FnTrixTri(Gs1Val,Gs2Val,Gs3Val,Gs1Pos,Gs2Pos,Gs3Pos,Ex1Val,Ex2Val,Ex1Pos,Ex2Pos )

%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
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
                        N2 = 6 - N1;
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
                                * nchoosek(N2,loopp)*(Ex1Val^loopp)*(Ex2Val^loopq)*sqrt(fp)*sqrt(fq);   
                        end
                    end
                end
            end
        end
    end
end

end

