clear all
close all
clc

% Display gaussians of BS used

% 1: He
% 2: H

r=0:0.01:4;

      ZETA1=2.0925;   % 2.0925
      ZETA2=1.24;     % 1.24

     COEF= [1.000000 0.000000 0.000000;
            0.678914 0.430129 0.000000;
            0.444635 0.535328 0.154329];
        
     EXPON=[0.270950 0.000000 0.000000;
            0.151623 0.851819 0.000000;
            0.109818 0.405771 2.227660];
 
sum_gauss1=zeros(3,length(r));
sum_gauss2=zeros(3,length(r));

for N=1:3
% SCALE THE EXPONENTS (A) OF PRIMITIVE GAUSSIANS
% INCLUDE NORMALIZATION IN CONTRACTION COEFFICIENTS (D)

sum_gauss1(N)=0;
sum_gauss2(N)=0;

figure(1)
subplot(3,2,1); title ('He, STO-1G')
subplot(3,2,2); title ('H,  STO-1G')
subplot(3,2,3); title ('He, STO-2G')
subplot(3,2,4); title ('H,  STO-2G')
subplot(3,2,5); title ('He, STO-3G')
subplot(3,2,6); title ('H,  STO-3G')

    for I=1:N
      % coefficients STO-NG for He
      A1(I)=EXPON(N,I)*(ZETA1^2);             %alpha
      D1(I)=COEF(N,I)*((2.0*A1(I)/pi)^0.75);  %unified coefficient: wf= SUM(D1*exp(-alpha*r^2))
      % coefficients STO-NG for H
      A2(I)=EXPON(N,I)*(ZETA2^2);
      D2(I)=COEF(N,I)*((2.0*A2(I)/pi)^0.75);
      
      gauss1(I,:)=D1(I)*exp(-A1(I)*r.^2);
      gauss2(I,:)=D2(I)*exp(-A2(I)*r.^2);
      
      subplot(3,2,2*N-1); hold on; plot(r,gauss1(I,:))
      subplot(3,2,2*N  ); hold on; plot(r,gauss2(I,:))
      
      sum_gauss1(N,:)=sum_gauss1(N,:)+gauss1(I,:);
      sum_gauss2(N,:)=sum_gauss2(N,:)+gauss2(I,:);
      
    end
    
      subplot(3,2,2*N-1); hold on; plot(r,sum_gauss1(N,:),'k','LineWidth',2)
      subplot(3,2,2*N  ); hold on; plot(r,sum_gauss2(N,:),'k','LineWidth',2)
   
end

%% Compare on the sampe plot

figure(2)
      subplot(1,2,1); hold on; plot(r,sum_gauss1); title ('He, STO-1,2,3G')
      subplot(1,2,2); hold on; plot(r,sum_gauss2); title ('H,  STO-1,2,3G')
