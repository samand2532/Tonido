function [ Created60Kets,Created42Kets,Created24Kets,Created06Kets ] =...
    FnBixTri(Gs1Val,Gs2Val,Gs1Pos,Gs2Pos,Ex1Val,Ex2Val,Ex3Val,Ex1Pos,Ex2Pos,Ex3Pos )
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
                for loopr = 0:6
                    if (loopi+loopj+loopp+loopq+loopr) == 6
                        
                        N1 = (loopi+loopj);
                        N2 = (loopp+loopq+loopr);
                        fN1 = factorial(N1);
                        fN2 = factorial(N2);
                        fi = factorial(loopi);
                        fj = factorial(loopj);
                        fp = factorial(loopp);
                        fq = factorial(loopq);
                        fr = factorial(loopr);
                        
                        if N1 == 6 && N2 == 0
                            aa = aa +1;
                            Created60Kets(aa,Gs1Pos+1) = loopi;
                            Created60Kets(aa,Gs2Pos+1) = loopj;
                            Created60Kets(aa,Ex1Pos+1) = loopp;
                            Created60Kets(aa,Ex2Pos+1) = loopq;
                            Created60Kets(aa,Ex3Pos+1) = loopr;
                            Created60Kets(aa,1) = (1/(sqrt(fN1*fN2))) * (fN2 / (fp*fq*fr))...
                                *(Gs1Val^loopi)*(Gs2Val^loopj)*sqrt(fi)*sqrt(fj)*nchoosek(N1,loopj)...
                                *(Ex1Val^loopp)*(Ex2Val^loopq)*(Ex3Val^loopr)*sqrt(fp)*sqrt(fq)*sqrt(fr);
                        elseif N1 == 4 && N2 == 2
                            bb = bb +1;
                            Created42Kets(bb,Gs1Pos+1) = loopi;
                            Created42Kets(bb,Gs2Pos+1) = loopj;
                            Created42Kets(bb,Ex1Pos+1) = loopp;
                            Created42Kets(bb,Ex2Pos+1) = loopq;
                            Created42Kets(bb,Ex3Pos+1) = loopr;
                            Created42Kets(bb,1) = (1/(sqrt(fN1*fN2))) * (fN2 / (fp*fq*fr))...
                                *(Gs1Val^loopi)*(Gs2Val^loopj)*sqrt(fi)*sqrt(fj)*nchoosek(N1,loopj)...
                                *(Ex1Val^loopp)*(Ex2Val^loopq)*(Ex3Val^loopr)*sqrt(fp)*sqrt(fq)*sqrt(fr);
                            elseif N1 == 2 && N2 == 4
                            cc = cc +1;
                            Created24Kets(cc,Gs1Pos+1) = loopi;
                            Created24Kets(cc,Gs2Pos+1) = loopj;
                            Created24Kets(cc,Ex1Pos+1) = loopp;
                            Created24Kets(cc,Ex2Pos+1) = loopq;
                            Created24Kets(cc,Ex3Pos+1) = loopr;
                            Created24Kets(cc,1) = (1/(sqrt(fN1*fN2))) * (fN2 / (fp*fq*fr))...
                                *(Gs1Val^loopi)*(Gs2Val^loopj)*sqrt(fi)*sqrt(fj)*nchoosek(N1,loopj)...
                                *(Ex1Val^loopp)*(Ex2Val^loopq)*(Ex3Val^loopr)*sqrt(fp)*sqrt(fq)*sqrt(fr);
                            elseif N1 == 0 && N2 == 6
                            dd = dd +1;
                            Created06Kets(dd,Gs1Pos+1) = loopi;
                            Created06Kets(dd,Gs2Pos+1) = loopj;
                            Created06Kets(dd,Ex1Pos+1) = loopp;
                            Created06Kets(dd,Ex2Pos+1) = loopq;
                            Created06Kets(dd,Ex3Pos+1) = loopr;
                            Created06Kets(dd,1) = (1/(sqrt(fN1*fN2))) * (fN2 / (fp*fq*fr))...
                                *(Gs1Val^loopi)*(Gs2Val^loopj)*sqrt(fi)*sqrt(fj)*nchoosek(N1,loopj)...
                                *(Ex1Val^loopp)*(Ex2Val^loopq)*(Ex3Val^loopr)*sqrt(fp)*sqrt(fq)*sqrt(fr);
                        end
                    end
                end
            end
        end
    end
end
                
                            
end

