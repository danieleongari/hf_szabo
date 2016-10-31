function INTGRL = INTGRL(IOP,N,R,ZETA1,ZETA2,ZA,ZB)
% C*********************************************************************
%       SUBROUTINE INTGRL(IOP,N,R,ZETA1,ZETA2,ZA,ZB)
% C
% C CALCULATES ALL THE BASIC INTEGRALS NEEDED FOR SCF CALCULATION
% C
% C*********************************************************************

global S12 T11 T12 T22 V11A V12A V22A V11B V12B V22B V1111 V2111 V2121 V2211 V2221 V2222

% THESE ARE THE CONTRACTION COEFFICIENTS AND EXPONENTS FOR
% A NORMALIZED SLATER ORBITAL WITH EXPONENT 1.0 IN TERMS OF
% NORMALIZED 1S PRIMITIVE GAUSSIANS
     
     COEF= [1.000000 0.000000 0.000000;
            0.678914 0.430129 0.000000;
            0.444635 0.535328 0.154329];
        
     EXPON=[0.270950 0.000000 0.000000;
            0.151623 0.851819 0.000000;
            0.109818 0.405771 2.227660];
     R2=R^2;
     
% SCALE THE EXPONENTS (A) OF PRIMITIVE GAUSSIANS
% INCLUDE NORMALIZATION IN CONTRACTION COEFFICIENTS (D)

for I=1:N
      % coefficients STO-NG for He
      A1(I)=EXPON(N,I)*(ZETA1^2);             %alpha
      D1(I)=COEF(N,I)*((2.0*A1(I)/pi)^0.75);  %unified coefficient
      % coefficients STO-NG for H
      A2(I)=EXPON(N,I)*(ZETA2^2);
      D2(I)=COEF(N,I)*((2.0*A2(I)/pi)^0.75);
end

% wf1_He = SUM(I=1,N) (D1(I)*exp(-A1(I)*r^2))
% wf2_H  = SUM(I=1,N) (D2(I)*exp(-A2(I)*r^2))
%
% D AND A ARE NOW THE CONTRACTION COEFFICIENTS AND EXPONENTS
% IN TERMS OF UNNORMALIZED PRIMITIVE GAUSSIANS

      S12=0.0;
      
      T11=0.0;
      T12=0.0;
      T22=0.0;
      
      V11A=0.0;
      V12A=0.0;
      V22A=0.0;
      V11B=0.0;
      V12B=0.0;
      V22B=0.0;
      
      V1111=0.0;
      V2111=0.0;
      V2121=0.0;
      V2211=0.0;
      V2221=0.0;
      V2222=0.0;
      
% CALCULATE ONE-ELECTRON INTEGRALS
% CENTER A IS FIRST ATOM, CETER B IS SECOND ATOM
% ORIGIN IS ON CENTER A
% V12A = OFF-DIAGONAL NUCLEAR ATTRACTION TO CENTER A, ETC.
for  I=1:N
    for J=1:N
     % RAP2 = SQUARED DISTANCE BETWEEN CENTER A AND CENTER P, ETC.
      RAP=A2(J)/(A1(I)+A2(J)) * R;
      RBP=A1(I)/(A1(I)+A2(J)) * R;
      RAP2=RAP^2;
      RBP2=RBP^2;
      
      S12=  S12   +S(A1(I),A2(J),R2)         *D1(I)*D2(J); % Overlapping wf1 and wf2
      T11=  T11   +T(A1(I),A1(J),0.0)        *D1(I)*D1(J); % Kinetic energy of isolated wf1 
      T12=  T12   +T(A1(I),A2(J),R2)         *D1(I)*D2(J); % Influenced K.en. between wf and the other nucleus
      T22=  T22   +T(A2(I),A2(J),0.0)        *D2(I)*D2(J); % Kinetic energy of isolated wf2
      V11A= V11A  +V(A1(I),A1(J),0.0,0.0,ZA) *D1(I)*D1(J); 
      V12A= V12A  +V(A1(I),A2(J),R2,RAP2,ZA) *D1(I)*D2(J);
      V22A= V22A  +V(A2(I),A2(J),0.0,R2,ZA)  *D2(I)*D2(J);
      V11B= V11B  +V(A1(I),A1(J),0.0,R2,ZB)  *D1(I)*D1(J);
      V12B= V12B  +V(A1(I),A2(J),R2,RBP2,ZB) *D1(I)*D2(J);
      V22B= V22B  +V(A2(I),A2(J),0.0,0.0,ZB) *D2(I)*D2(J);
      
    end
end

% CALCULATE TWO-ELECTRON INTEGRALS
for I=1:N
    for J=1:N
        for K=1:N
            for L=1:N
      RAP=A2(I)*R/(A2(I)+A1(J));
      RBP=R-RAP;
      RAQ=A2(K)*R/(A2(K)+A1(L));
      RBQ=R-RAQ;
      RPQ=RAP-RAQ;
      RAP2=RAP*RAP;
      RBP2=RBP*RBP;
      RAQ2=RAQ*RAQ;
      RBQ2=RBQ*RBQ;
      RPQ2=RPQ*RPQ;
      
      V1111= V1111 +TWOE(A1(I),A1(J),A1(K),A1(L),0.0,0.0, 0.0) *D1(I)*D1(J)*D1(K)*D1(L);
      V2111= V2111 +TWOE(A2(I),A1(J),A1(K),A1(L), R2,0.0,RAP2) *D2(I)*D1(J)*D1(K)*D1(L);
      V2121= V2121 +TWOE(A2(I),A1(J),A2(K),A1(L), R2, R2,RPQ2) *D2(I)*D1(J)*D2(K)*D1(L);
      V2211= V2211 +TWOE(A2(I),A2(J),A1(K),A1(L),0.0,0.0,  R2) *D2(I)*D2(J)*D1(K)*D1(L);
      V2221= V2221 +TWOE(A2(I),A2(J),A2(K),A1(L),0.0, R2,RBQ2) *D2(I)*D2(J)*D2(K)*D1(L);
      V2222= V2222 +TWOE(A2(I),A2(J),A2(K),A2(L),0.0,0.0, 0.0) *D2(I)*D2(J)*D2(K)*D2(L);
      
            end
        end
    end
end

    if IOP~=0
        disp(['R = ',num2str(R)])
        disp(['ZETA1 = ',num2str(ZETA1)])
        disp(['ZETA2 = ',num2str(ZETA2)])
        disp(['S12 = ',num2str(S12)])
        disp(['T11 = ',num2str(T11)])
        disp(['T12 = ',num2str(T12)])
        disp(['T22 = ',num2str(T22)])
        disp(['V11A = ',num2str(V11A)])
        disp(['V12A = ',num2str(V12A)])
        disp(['V22A = ',num2str(V22A)])
        disp(['V11B = ',num2str(V11B)])
        disp(['V12B = ',num2str(V12B)])
        disp(['V22B = ',num2str(V22B)])
        disp(['V1111 = ',num2str(V1111)])
        disp(['V2111 = ',num2str(V2111)])
        disp(['V2121 = ',num2str(V2121)])
        disp(['V2211 = ',num2str(V2211)])
        disp(['V2221 = ',num2str(V2221)])
        disp(['V2222 = ',num2str(V2222)])
    end


end

