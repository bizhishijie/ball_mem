% %%
% % 定义常数 国际单位
% a = 0.047;% 鼓的半径
% ca = 346;% 空气中声速
% cs = sqrt(7/0.036);% 橡胶中波的传播速度
% L = 0.117;% 鼓的深度
% gamma = 1.4; % 空气绝热系数
% rho = 1.184;% 空气密度
% % d = 3e-4; % 气球厚度
% sigma = 0.036;% 气球面密度
% pa = 1e5;% 大气压强
% 
% cs_ca=cs/ca;
% L_a=L/a;
% sigma_a_rho=sigma/a/rho;
% cs_ca=0.3;sigma_a_rho=3.3;L_a=1.5;
%%
a =0.047;% 鼓的半径
L = 0.1;% 鼓的深度
T=25; %张力19.9
gamma = 1.4; % 空气绝热系数
rho = 1.184;% 空气密度
% d = 3e-4; % 气球厚度
ca = 345;% 空气中声速
sigma = 0.095;% 气球面密度
pa = 1e5;% 大气压强
cs = sqrt(T/sigma);% 橡胶中波的传播速度
cs_ca=cs/ca;L_a=L/a;sigma_a_rho=sigma/a/rho;
%%
load('rootBessel.mat')
load('rootBesselDiff.mat')
for n=0:10
    disp(n)
    OmegaResonance(n,10, cs_ca,L_a, sigma_a_rho ,rootBessel,rootBesselDiff);
end
