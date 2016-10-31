function V = V(A,B,RAB2,RCP2,ZC)
% C*********************************************************************
%       FUNCTION V(A,B,RAB2,RCP2,ZC)
% C
% C CALCULATES UN-NORMALIZED NUCLEAR ATTRACTION INTEGRALS
% C
% C*********************************************************************

V=2*pi/(A+B)* F0((A+B)*RCP2)*exp(-A*B*RAB2/(A+B));
V=-V*ZC;

end

