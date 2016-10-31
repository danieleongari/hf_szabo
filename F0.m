function F0 = F0(ARG)
% C*********************************************************************
%       FUNCTION F0(ARG)
% C
% C CALCULATES THE F FUNCTION
% C FO ONLY (S-TYPE ORBITALS)
% C
% C*********************************************************************

if ARG>=1.0D-6  %F0 IN TERMS OF THE ERROR FUNCTION
    
    F0=sqrt(pi/ARG)*erf(sqrt(ARG))/2;
      
else %ASYMPTOTIC VALUE FOR SMALL ARGUMENTS
    
    F0=1-ARG/3;
    
end

end

