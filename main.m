clear all
close all
clc

      IOP=2;                    
      N=3;     
      R=1.4632;       % 1.4632
      ZETA1=2.0925;   % 2.0925
      ZETA2=1.24;     % 1.24
      ZA=2;           % 2
      ZB=1;           % 1 
      HFCALC(IOP,N,R,ZETA1,ZETA2,ZA,ZB); %explanation of all parameters     
      
     
      
%% E(R) (alternative to previous one: put it in ctrl+R )

% global etot eelec enucl ri
% 
%       IOP=1;%only converged results                
%       N=3;     
%       ZETA1=2.0925;   % 2.0925
%       ZETA2=1.24;     % 1.24
%       ZA=2;           % 2
%       ZB=1;           % 1 
%      
%       r_scan_min=0.5;
%       r_scan_max=2.5;
%       r_scan_step=0.1;
%       
%       ri=0;
%       for r=r_scan_min:r_scan_step:r_scan_max
%         ri=ri+1;
%         rplot(ri)=r;
%         HFCALC(IOP,N,r,ZETA1,ZETA2,ZA,ZB); %explanation of all parameters 
%       end
%       
%       plot(rplot,etot,'b')
%       hold on
%       plot(rplot,eelec,'r')
%       plot(rplot,enucl,'g')
%       legend('E_{TOT}','E_{ELECT}','E_{NUCL}','Location','NorthEastOutside')
      
