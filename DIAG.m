function [C,E] = DIAG(F,C,E)
% C^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^*
%       SUBROUTINE DIAG(F,C,E)
% C
% C DIAGONALIZES F TO GIVE EIGENVECTORS IN C AND EIGENVALUES IN E
% C THETA IS THE ANGLE DESCRIBING SOLUTION
% C
% C^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^*


if abs(F(1,1)-F(2,2)) < 1.0D-20
%C HERE IS SYMMETRY DETERMINED SOLUTION (HOMONUCLEAR DIATOMIC)
      THETA=pi/4;
else
%C SOLUTION FOR HETERONUCLEAR DIATOMIC
      THETA=0.5*atan(2*F(1,2)/(F(1,1)-F(2,2)));
end

      C(1,1)=cos(THETA);
      C(2,1)=sin(THETA);
      C(1,2)=sin(THETA);
      C(2,2)=-cos(THETA);
      E(1,1)=F(1,1)*cos(THETA)^2+F(2,2)*sin(THETA)^2+F(1,2)*sin(2*THETA);
      E(2,2)=F(2,2)*cos(THETA)^2+F(1,1)*sin(THETA)^2-F(1,2)*sin(2*THETA);
      E(2,1)=0;
      E(1,2)=0;
      
%C ORDER EIGENVALUES AND EIGENVECTORS
if E(2,2)<E(1,1)
      TEMP=E(2,2);
      E(2,2)=E(1,1);
      E(1,1)=TEMP;
      TEMP=C(1,2);
      C(1,2)=C(1,1);
      C(1,1)=TEMP;
      TEMP=C(2,2);
      C(2,2)=C(2,1);
      C(2,1)=TEMP;
end

end

