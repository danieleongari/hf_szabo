close all
clear all
clc

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Program that compares the energy for the hydrogen atom between:
% - the exact solution
% - STO-1G, STO-2G and STO-3G
% printing the norm (I), the kinetic (K), the potential (V) and the total (E) energy.
% BC: requires symbolic package (syms) of MatLab and therefore is not compatible with Octave.
%     Therefore, you can find the screenshot of the solution as hydrogen_1S.png to compare.
% (Daniele, 16 November 2016)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% 1S orbital
n=1;
l=0;

% http://hyperphysics.phy-astr.gsu.edu/hbase/quantum/hydwf.html#c3
hbar=1; % plank_constant/(2*pi)
mass=1; % electron_mass
e=1;    % |electron_charge|

a0=hbar^2/mass/e^2; % Bohr_radius (0.529 angstrom)

% Exact solution (spherical coordinates)-----------------------------------
syms phi theta r
F=1/sqrt(2*pi);            % F(phi)
P=1/sqrt(2);               % P(theta)
R=2/a0^(3/2) * exp(-r/a0); % R(r)

wf=F*P*R; %wave_function(fi,theta,r)
radial_prob=int(int(wf^2*r*sin(theta),theta,0,pi),phi,0,2*pi);

I=double(int(int(int(wf^2*r^2*sin(theta),r,0,20),theta,0,pi),phi,0,2*pi));

Kwf=-0.5*r^(-2)*diff(r^2*diff(wf,r),r);
K=double(int(int(int(wf*Kwf*r^2*sin(theta),r,0,20),theta,0,pi),phi,0,2*pi));
V=double(int(int(int(wf*(-1/r)*wf*r^2*sin(theta),r,0,20),theta,0,pi),phi,0,2*pi));
E=K+V;

%plot
fprintf('Exact solution: I= %f, K= %f, V= %f, E=K+V= %f \n',I,K,V,E)
subplot(2,2,1)
ezplot(wf,[0,4,0,1])
hold on
ezplot(radial_prob,[0,4,0,1])
title(sprintf('Exact solution (E=%f)',E))
legend('wf(r)','Radial prob.(r)')

%Plot best solution for one gaussian---------------------------------------
% see Szabo Ostlund ex 1.19

alpha=8/(9*pi); %strangely STO-1G in the hf_szabo program is different
N=(2*alpha/pi)^(3/4);
gauss1=N*exp(-alpha*r^2);
wf1=gauss1;
I=double(int(int(int(wf1^2*r^2*sin(theta),r,0,20),theta,0,pi),phi,0,2*pi));

Kwf=-0.5*r^(-2)*diff(r^2*diff(wf1,r),r);
K=double(int(int(int(wf1*Kwf*r^2*sin(theta),r,0,20),theta,0,pi),phi,0,2*pi));
V=double(int(int(int(wf1*(-1/r)*wf1*r^2*sin(theta),r,0,20),theta,0,pi),phi,0,2*pi));
E=K+V;

%plot
fprintf('Single gauss: I= %f, K= %f, V= %f, E=K+V= %f \n',I,K,V,E)
subplot(2,2,2)
ezplot(wf,[0,4,0,1])
hold on
ezplot(wf1,[0,4,0,1])
title(sprintf('Single G (E=%f)',E))
legend('wf(r)','gauss1(r)')

%STO-2G (Szabo Ostlund = Hehre,Pople 1969 = BSexchange)--------------------
zeta=1.24; %for H
alpha=[0.151623,0.851819]*zeta^2;
d=[0.678914,0.430129];
sto2g=0;
for i=1:2
    N=(2*alpha(i)/pi)^(3/4);
    sto2g=sto2g+d(i)*N*exp(-alpha(i)*r^2);
end
wf1=sto2g;
I=double(int(int(int(wf1^2*r^2*sin(theta),r,0,20),theta,0,pi),phi,0,2*pi));

Kwf=-0.5*r^(-2)*diff(r^2*diff(wf1,r),r);
K=double(int(int(int(wf1*Kwf*r^2*sin(theta),r,0,20),theta,0,pi),phi,0,2*pi));
V=double(int(int(int(wf1*(-1/r)*wf1*r^2*sin(theta),r,0,20),theta,0,pi),phi,0,2*pi));
E=K+V;

%plot
fprintf('STO-2G: I= %f, K= %f, V= %f, E=K+V= %f \n',I,K,V,E)
subplot(2,2,3)
ezplot(wf,[0,4,0,1])
hold on
ezplot(wf1,[0,4,0,1])
title(sprintf('STO-2G (E=%f)',E))
legend('wf(r)','STO-2G(r)')

%STO-3G (Szabo Ostlund = Hehre,Pople 1969 = BSexchange)--------------------
zeta=1.24; %for H
alpha=[0.109818,0.405771,2.227660]*zeta^2;
d=[0.444635,0.535328,0.154329];
sto3g=0;
for i=1:3
    N=(2*alpha(i)/pi)^(3/4);
    sto3g=sto3g+d(i)*N*exp(-alpha(i)*r^2);
end
wf1=sto3g;
I=double(int(int(int(wf1^2*r^2*sin(theta),r,0,20),theta,0,pi),phi,0,2*pi));

Kwf=-0.5*r^(-2)*diff(r^2*diff(wf1,r),r);
K=double(int(int(int(wf1*Kwf*r^2*sin(theta),r,0,20),theta,0,pi),phi,0,2*pi));
V=double(int(int(int(wf1*(-1/r)*wf1*r^2*sin(theta),r,0,20),theta,0,pi),phi,0,2*pi));
E=K+V;

%plot
fprintf('STO-2G: I= %f, K= %f, V= %f, E=K+V= %f \n',I,K,V,E)
subplot(2,2,4)
ezplot(wf,[0,4,0,1])
hold on
ezplot(wf1,[0,4,0,1])
title(sprintf('STO-3G (E=%f)',E))
legend('wf(r)','STO-3G(r)')
