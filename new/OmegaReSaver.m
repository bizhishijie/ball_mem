%%
% 定义常数 国际单位
% a = 1;% 鼓的半径
% ca = 340;% 空气中声速
% cs = 0.6*ca;% 橡胶中波的传播速度
% L = 5*a;% 鼓的深度
% gamma = 1.4; % 空气绝热系数
% rho = 1.5e3;%橡胶密度
% d = 1e-4; % 气球厚度
% sigma = 3.3*a*rho;% 气球面密度
% pa = 1e5;% 大气压强
%%
load('rootBessel.mat')
load('rootBesselDiff.mat')
for n=0:50
    disp(n)
    OmegaResonanceNew(n,50, 0.6,5, 3.3 ,rootBessel,rootBesselDiff)
end