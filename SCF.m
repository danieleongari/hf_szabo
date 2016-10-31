function SCF = SCF(IOP,N,R,ZETA1,ZETA2,ZA,ZB)
% C*********************************************************************
%       SUBROUTINE SCF(IOP,N,R,ZETA1,ZETA2,ZA,ZB)
% C
% C PERFORMS THE SCF ITERATIONS
% C
% C*********************************************************************

global s X XT H F G C FPRIME CPRIME P OLDP TT E
global etot eelec enucl ri

%%  parameters
%C CONVERGENCE CRITERION FOR DENSITY MATRIX
CRIT=1.0D-4;

%C MAXIMUM NUMBER OF ITERATIONS
MAXIT=25;

%%
%C ITERATION NUMBER
    ITER=0;
    
%C USE CORE-HAMILTONIAN FOR INITIAL GUESS AT F, I.E. (P=0)
      P=zeros(2,2);
if IOP==2
      P
end

%C START OF ITERATION LOOP
while ITER<MAXIT
    
ITER=ITER+1;
    
if IOP==2 
    disp(['-------------- START OF ITERATION NUMBER = ',num2str(ITER)])
end

%C FORM TWO-ELECTRON PART OF FOCK MATRIX FROM P
%----------------CALL FORMG------------------
% C CALCULATES THE G MATRIX FROM THE DENSITY MATRIX
% C AND TWO-ELECTRON INTEGRALS      
G=zeros(2,2);

for I=1:2
for J=1:2
for K=1:2
for L=1:2      
      G(I,J)=G(I,J)+P(K,L)*(TT(I,J,K,L)-0.5*TT(I,L,K,J));
end
end
end
end  
      
if IOP==2
	G
end

%C ADD CORE HAMILTONIAN TO GET FOCK MATRIX 
    F = H+G; %(2,2)matrix

%C CALCULATE ELECTRONIC ENERGY
    EN=0.0;
for I=1:2
for J=1:2
    EN=EN+0.5*P(I,J)*(H(I,J)+F(I,J));
end
end

if IOP==2 
	F
   	disp(['ELECTRONIC ENERGY = ',num2str(EN)])
end

%C TRANSFORM FOCK MATRIX USING G FOR TEMPORARY STORAGE
      G=F*X; %(2,2) matrix
      FPRIME=XT*G; %(2,2) matrix
      
%C DIAGONALIZE TRANSFORMED FOCK MATRIX
      [CPRIME,E]=DIAG(FPRIME,CPRIME,E);
      
%C TRANSFORM EIGENVECTORS TO GET MATRIX C
      C=X*CPRIME; %(2,2) matrix
      
%C FORM NEW DENSITY MATRIX
      OLDP=P; %C SAVE PRESENT DENSITY MATRIX BEFORE CREATING NEW ONE
      P=zeros(2,2);
      
 for I=1:2
 for J=1:2
 for K=1
      P(I,J)=P(I,J)+2.0D0*C(I,K)*C(J,K);
 end
 end
 end
  
if IOP==2
    FPRIME
    CPRIME
    E
    C
    P
end

%C CALCULATE DELTA
      DELTA=0.0;
for I=1:2
for J=1:2
    DELTA=DELTA+(P(I,J)-OLDP(I,J))^2;
end
end
      DELTA=sqrt(DELTA/4);
if IOP~=0

	disp(['DELTA(CONVERGENCE OF DENSITY MATRIX) = ',num2str(DELTA)])
    disp ' '
end

%C CHECK FOR CONVERGENCE
if DELTA<CRIT

    %C CALCULATION CONVERGED IF IT GOT HERE
    %C ADD NUCLEAR REPULSION TO GET TOTAL ENERGY

    ENT=EN+ZA*ZB/R;
    
                     if exist('ri') % scan E(R)
                       etot(ri)=ENT;
                       eelec(ri)=EN;
                       enucl(ri)=ZA*ZB/R;
                     end
    
    if IOP~=0
        disp 'CALCULATION CONVERGED'
        disp ''
        disp(['ELECTRONIC ENERGY = ',num2str(EN)])
        disp(['TOTAL ENERGY = ',num2str(ENT)])
    end

    if IOP==1
    % C PRINT OUT THE FINAL RESULTS IF
    % C HAVE NOT DONE SO ALREADY
        G
        F
        E
        C
        P
    end

    %C PS MATRIX HAS MULLIKEN POPULATIONS
      OLDP=P*s; %(2,2) matrix

    if IOP~=0
        OLDP
    end
    
    return %break the SCF function
end

end

%C ITER=MAXIT
%C SOMETHING WRONG HERE
disp 'NO CONVERGENCE IN SCF AFTER MAXITER'

end

