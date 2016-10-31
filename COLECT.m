function COLECT = COLECT(IOP,N,R,ZETA1,ZETA2,ZA,ZB)
% C*********************************************************************
%       SUBROUTINE COLECT(IOP,N,R,ZETA1,ZETA2,ZA,ZB)
% C
% C THIS TAKES THE BASIC INTEGRALS FROM COMMON AND ASSEMBLES THE
% C RELEVENT MATRICES, THAT IS s,H,X,XT, AND TWO-ELECTRON INTEGRALS
% C
% C*********************************************************************

global S12 T11 T12 T22 V11A V12A V22A V11B V12B V22B V1111 V2111 V2121 V2211 V2221 V2222
global s X XT H F G C FPRIME CPRIME P OLDP TT E

%declaration of all matrix used in COLECT
s=      zeros(2,2);
X=      zeros(2,2);
XT=     zeros(2,2);
H=      zeros(2,2);
F=      zeros(2,2);
G=      zeros(2,2);
C=      zeros(2,2);
FPRIME= zeros(2,2);
CPRIME= zeros(2,2);
P=      zeros(2,2);
OLDP=   zeros(2,2);
TT=     zeros(2,2,2,2);%TT is a 4dimensional matrix!
E=      zeros(2,2);

%C FORM CORE HAMILTONIAN
      H(1,1)=T11+V11A+V11B;
      H(1,2)=T12+V12A+V12B;
      H(2,1)=H(1,2);
      H(2,2)=T22+V22A+V22B;
%C FORM OVERLAP MATRIX
      s(1,1)=1;
      s(1,2)=S12;
      s(2,1)=S12;
      s(2,2)=1;
%C USE CANONICAL ORTHOGONALIZATION
      X(1,1)=1/sqrt(2*(1+S12));
      X(2,1)=X(1,1);
      X(1,2)=1/sqrt(2*(1-S12));
      X(2,2)=-X(1,2);
%C TRANSPOSE OF TRANSFORMATION MATRIX
      XT(1,1)=X(1,1);
      XT(1,2)=X(2,1);
      XT(2,1)=X(1,2);
      XT(2,2)=X(2,2);
%C MATRIX OF TWO-ELECTRON INTEGRALS
      TT(1,1,1,1)=V1111; 
      TT(2,1,1,1)=V2111;
      TT(1,2,1,1)=V2111;
      TT(1,1,2,1)=V2111;
      TT(1,1,1,2)=V2111;
      TT(2,1,2,1)=V2121;
      TT(1,2,2,1)=V2121;
      TT(2,1,1,2)=V2121;
      TT(1,2,1,2)=V2121;
      TT(2,2,1,1)=V2211;
      TT(1,1,2,2)=V2211;
      TT(2,2,2,1)=V2221;
      TT(2,2,1,2)=V2221;
      TT(2,1,2,2)=V2221;
      TT(1,2,2,2)=V2221;
      TT(2,2,2,2)=V2222;
      
if IOP~=0
      s
      X
      H
      
    disp ' '
    disp 'TT='
    disp ' '
    for I=1:2
    for J=1:2
    for K=1:2
    for L=1:2
    disp(['( ',num2str(I),' ',num2str(J),' ',num2str(K),' ',num2str(L),' ) ',num2str(TT(I,J,K,L))])
    end
    end
    end
    end
end

end

