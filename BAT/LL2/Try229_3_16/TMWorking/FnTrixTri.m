function [ Created60Kets,Created42Kets,Created24Kets,Created06Kets ] =...
    FnTrixTri(Gs1Val,Gs2Val,Gs3Val,Gs1Pos,Gs2Pos,Gs3Pos,Ex1Val,Ex2Val,Ex3Val,Ex1Pos,Ex2Pos,Ex3Pos )
%The TrixTri expansion, the Basis and Conts Created.


Created60Kets = zeros(210,23);
Created42Kets = zeros(210,23);
Created24Kets = zeros(210,23);
Created06Kets = zeros(210,23);
aa=0;bb=0;cc=0;dd=0;
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

end

