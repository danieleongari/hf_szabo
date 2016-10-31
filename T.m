function T = T(A,B,RAB2)
% C*********************************************************************
%       FUNCTION T(A,B,RAB2)
% C
% C CALCULATES KINETIC ENERGY INTEGRALS FOR UN-NORMALIZED PRIMITIVES
% C
% C*********************************************************************

T=A*B/(A+B)*(3-2*A*B*RAB2/(A+B))*(pi/(A+B))^1.5*exp(-A*B*RAB2/(A+B));

end

