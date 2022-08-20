%%
% % 定义常数 国际单位
a=0.047; %半径
L = 0.1;% 鼓的深度
gamma = 1.4; % 空气绝热系数
% T=20; % 张力19.9
rho = 1.184;% 空气密度
ca = 345;% 空气中声速
sigma = 0.120;% 气球面密度
pa = 1e5;% 大气压强
cs = sqrt(T/sigma);% 橡胶中波的传播速度
cs_ca=cs/ca;L_a=L/a;sigma_a_rho=sigma/a/rho;
load('rootBessel.mat')
load('rootBesselDiff.mat')
%%
% L_a=2;sigma_a_rho=3.3;cs_ca=0.653185;
%%
parfor n=0:20
    disp(n)
    OmegaResonance(n,20, cs_ca,L_a, sigma_a_rho,rootBessel,rootBesselDiff);
end