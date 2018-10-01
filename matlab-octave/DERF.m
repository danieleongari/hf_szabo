function DERF = DERF(ARG)
% C*********************************************************************
%       FUNCTION DERFOTHER(ARG)
% C
% C CALCULATES THE ERROR FUNCTION ACCORDING TO A RATIONAL
% C APPROXIMATION FROM M. ARBRAMOWITZ AND I.A. STEGUN,
% C HANDBOOK OF MATHEMATICAL FUNCTIONS, DOVER.
% C ABSOLUTE ERROR IS LESS THAN 1.5*10**(-7)
% C CAN BE REPLACED BY A BUILT-IN FUNCTION ON SOME MACHINES
% C
% C Daniele: unused, I'm using directly MatLab erf(x) function
% C
% C*********************************************************************

      P=0.3275911D0;
      A=[0.254829592D0,-0.284496736D0,1.421413741D0,-1.453152027D0,1.061405429D0];
      t=1.0D0/(1.0D0+P*ARG); %t instead of T to avoid confusion with function T
      TN=t;
      
      POLY=A(1)*TN;
 for I=2:5
      TN=TN*t;
      POLY=POLY+A(I)*TN;
 end
 
      DERF=1.0D0-POLY*exp(-ARG*ARG);
      
end

