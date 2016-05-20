function [ Created60Kets,Created42Kets,Created24Kets,Created06Kets ] =...
    FnBixBi(Gs1Val,Gs2Val,Gs1Pos,Gs2Pos,Ex1Val,Ex2Val,Ex1Pos,Ex2Pos )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
Created60Kets = zeros(210,23);
Created42Kets = zeros(210,23);
Created24Kets = zeros(210,23);
Created06Kets = zeros(210,23);
aa=0;bb=0;cc=0;dd=0;ee=0;
for loopi = 0:6
    for loopj = 0:6
        for loopp = 0:6
            for loopq = 0:6
                
                    if (loopi+loopj+loopp+loopq) == 6
                        
                        N1 = (loopi+loopj);
                        N2 = (loopp+loopq);
                        fN1 = factorial(N1);
                        fN2 = factorial(N2);
                        fi = factorial(loopi);
                        fj = factorial(loopj);
                        fp = factorial(loopp);
                        fq = factorial(loopq);
                        
                        
                        if N1 == 6 && N2 == 0
                            aa = aa +1;
                            Created60Kets(aa,Gs1Pos+1) = loopi;
                            Created60Kets(aa,Gs2Pos+1) = loopj;
                            Created60Kets(aa,Ex1Pos+1) = loopp;
                            Created60Kets(aa,Ex2Pos+1) = loopq;
                            
                            Created60Kets(aa,1) = (1/(sqrt(fN1*fN2)))...
                                *(Gs1Val^loopi)*(Gs2Val^loopj)*sqrt(fi)*sqrt(fj)*nchoosek(N1,loopj)...
                                *(Ex1Val^loopp)*(Ex2Val^loopq)*sqrt(fp)*sqrt(fq)*nchoosek(N2,loopq);
                        elseif N1 == 4 && N2 == 2
                            bb = bb +1;
                            Created42Kets(bb,Gs1Pos+1) = loopi;
                            Created42Kets(bb,Gs2Pos+1) = loopj;
                            Created42Kets(bb,Ex1Pos+1) = loopp;
                            Created42Kets(bb,Ex2Pos+1) = loopq;
                            Created42Kets(bb,1) = (1/(sqrt(fN1*fN2)))...
                                *(Gs1Val^loopi)*(Gs2Val^loopj)*sqrt(fi)*sqrt(fj)*nchoosek(N1,loopj)...
                                *(Ex1Val^loopp)*(Ex2Val^loopq)*sqrt(fp)*sqrt(fq)*nchoosek(N2,loopq);
                            elseif N1 == 2 && N2 == 4
                            cc = cc +1;
                            Created24Kets(cc,Gs1Pos+1) = loopi;
                            Created24Kets(cc,Gs2Pos+1) = loopj;
                            Created24Kets(cc,Ex1Pos+1) = loopp;
                            Created24Kets(cc,Ex2Pos+1) = loopq;
                            
                            Created24Kets(cc,1) = (1/(sqrt(fN1*fN2)))...
                                *(Gs1Val^loopi)*(Gs2Val^loopj)*sqrt(fi)*sqrt(fj)*nchoosek(N1,loopj)...
                                *(Ex1Val^loopp)*(Ex2Val^loopq)*sqrt(fp)*sqrt(fq)*nchoosek(N2,loopq);
                            elseif N1 == 0 && N2 == 6
                            dd = dd +1;
                            Created06Kets(dd,Gs1Pos+1) = loopi;
                            Created06Kets(dd,Gs2Pos+1) = loopj;
                            Created06Kets(dd,Ex1Pos+1) = loopp;
                            Created06Kets(dd,Ex2Pos+1) = loopq;
                            
                            Created06Kets(dd,1) = (1/(sqrt(fN1*fN2)))...
                                *(Gs1Val^loopi)*(Gs2Val^loopj)*sqrt(fi)*sqrt(fj)*nchoosek(N1,loopj)...
                                *(Ex1Val^loopp)*(Ex2Val^loopq)*sqrt(fp)*sqrt(fq)*nchoosek(N2,loopq);
                        end
                    end
                end
            end
        end
    end
end



