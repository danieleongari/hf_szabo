function S = S(A,B,RAB2)
% C*********************************************************************
%       FUNCTION S(A,B,RAB2)
% C
% C CALCULATES OVERLAPS FOR UN-NORMALIZED PRIMITIVES
% C
% C*********************************************************************

S=(pi/(A+B))^1.5*exp(-A*B*RAB2/(A+B));

end

