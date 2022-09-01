% a=0.05;
% L=0.1;
% rho_mem=0.095;
% rho_air=1.184;
% a0=0.1*a;
% T=41.9;
% cs=sqrt(T/rho_mem);
% ca=345;
% 
% sigma = 0.120;

% cs_ca=cs/ca;L_a=L/a;sigma_a_rho=sigma/a/rho_air;

L_a=5;
sigma_a_rho=3.3;
cs_ca=0.6;
a=1;
ca=345;cs=cs_ca*ca;

load('rootBessel.mat')
load('rootBesselDiff.mat')

parfor n=0:10
    disp(n)
    OmegaResonance(n,10, cs_ca,L_a, sigma_a_rho,rootBessel,rootBesselDiff);
end
OmegaReLoader
load('OmegaRe.mat');

rootBesselDiff=rootBesselDiff(1,:);
% rootBessel=rootBessel(1,:);
Omega=Omega(1,:);
f_cav_close=sqrt(rootBesselDiff.^2+pi^2)/L_a;
f_cav_open=sqrt(rootBesselDiff.^2+(pi/2)^2)/L_a;
f_mem=rootBessel*cs_ca;
f_cou=Omega;

f_cav_close=sort(f_cav_close(:));
f_cav_open=sort(f_cav_open(:));
f_mem=sort(f_mem(:));
f_cou=sort(f_cou(:));

clf;hold on

plot(f_cav_close,4*ones(size(f_cav_close)),'k*')
plot(f_cav_open,3*ones(size(f_cav_open)),'g*')
plot(f_mem,2*ones(size(f_mem)),'r*')
plot(f_cou,ones(size(f_cou)),'b*')

axis([0.3 10 0 5]);
grid on
set(gca,'xscale','log')