clear all;
  InitialBras = [5 0 0 0 0 0 1 ;...
                 4 1 0 0 0 1 0 ;...
                 4 0 1 0 1 0 0 ;...
                 3 2 0 0 1 0 0 ;...
                 4 0 0 2 0 0 0 ;...
                 3 1 1 1 0 0 0 ;...
                 2 3 0 1 0 0 0 ;...
                 3 0 3 0 0 0 0 ;...
                 2 2 2 0 0 0 0 ;...
                 1 4 1 0 0 0 0 ;...
                 0 6 0 0 0 0 0 ];



  InitialKets = [5 4 4 3 4 3 2 3 2 1 0;...
                 0 1 0 2 0 1 3 0 2 4 6;...
                 0 0 1 0 0 1 0 3 2 1 0;...
                 0 0 0 0 2 1 1 0 0 0 0;...
                 0 0 1 1 0 0 0 0 0 0 0;...
                 0 1 0 0 0 0 0 0 0 0 0;...
                 1 0 0 0 0 0 0 0 0 0 0];
 CreaVect = zeros(7,1);
 AnnaVect = zeros(7,1);  
 P1CreaVect = CreaVect; P2CreaVect = CreaVect;
 P1AnnaVect = AnnaVect; P2AnnaVect = AnnaVect;
             
             % m cannot go < 1 or >7
             %REMEMBER THAT HAVE TO START AT 1!!!
             % p1 = root(m(m-1))*adagm*acrea(m-2)
             % p2 = root((m+1)(m+2))*adagm * acrea (m+2)

  bra = 1;
  ket = 1;

  InitialBra = [1 1 1 1 1 1 1]'%InitialBras(bra,:)';
  InitialKet = [1 1 1 1 1 1 1]'%InitialKets(:,ket);
  P1=0;
  P2=0;
  

  m = 0 % from 1 to 7 - FOR MATLAB - ann / crea changed 
  if (m-1)>1
      P1CreaVect(m-1,1) = 1;
  else
      P1CreaVect = InitialKet.*9; %%%% .* 9 to ensure non-orthogonal
  end
  
  if (m-3)>0
      P1AnnaVect(m-3,1) = -1;
      else
      P1AnnaVect = InitialKet.*9; %%%% .* 9 to ensure non-orthogonal
  end
  if (m-1)>1
      P2CreaVect(m-1,1) = 1;
      else
      P2CreaVect = InitialKet.*9; %%%% .* 9 to ensure non-orthogonal
  end
  P2AnnaVect(m+1,1) = -1;

  %%%%%P1
  if (m-3) >=0
      P1=sqrt((m-1)*(m-2))*sqrt(InitialBra((m-1),1))*sqrt(InitialBra((m-3),1));
      P1trans = InitialKet + P1CreaVect + P1AnnaVect;
  else
      P1=0;
      P1trans = InitialKet .*9;
  end
  
  
  
      
  
  
  
  
  P1
  %P2
  P1trans
  %P2trans
  




