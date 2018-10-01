function TWOE = TWOE(A,B,C,D,RAB2,RCD2,RPQ2)
% C**********************************************************************
%       FUNCTION TWOE(A,B,C,D,RAB2,RCD2,RPQ2)
% C
% C CALCULATES TWO-ELECTRON INTEGRALS FOR UN-NORMALIZED PRIMITIVES
% C A,B,C,D ARE THE EXPONENTS ALPHA, BETA, ETC.
% C RAB2 EQUALS SQUARED DISTANCE BETWEEN CENTER A AND CENTER B, ETC.
% C**********************************************************************

TWOE=2*(pi^2.5)/((A+B)*(C+D)*sqrt(A+B+C+D))*F0((A+B)*(C+D)*RPQ2/(A+B+C+D))*exp(-A*B*RAB2/(A+B)-C*D*RCD2/(C+D));

end

