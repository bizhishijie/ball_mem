%%
% 定义常数 国际单位
a=0.2;% 鼓的半径
ca=340;% 空气中声速
cs=1800;% 橡胶中声速
L=0.3;% 鼓的深度
    gamma=1.4; % 空气绝热系数
ro=1.5e3;%橡胶密度
d=1e-4; % 气球厚度
sigma=ro*d;% 气球面密度
%%
% 计算单独膜的频率
% r,m,n
r0=1;%(15) r=1
load('rootBessel.mat','rootBessel');%读取贝塞尔方程的近似根
rDiff=100;%将r划分的网格数
r=(0:r0/(rDiff-1):r0)';
m=1:10;
n=0:10;%划分网格

mu=permute(rootBessel,[3,2,1]);
mu=repmat(mu,length(r),1,1);% 重整mu

sqrt2=sqrt(2);
phi=zeros(size(mu));
mur=mu.*r;
for ni=n+1
    phi(:,:,ni)=sqrt2*besselj(ni-1,mur(:,:,ni))./besselj(ni-1,mu(:,:,ni));
end
%%
% 画关于r的函数(5)
% plot(r,phi(:,3,1))% 数值特别大，在r=1处为0
% 对其归一化
sum1=sum(phi.*phi.*r,1)*r0/rDiff;% (7)
phi=phi./sqrt(repmat(sum1,length(r),1,1));
%%
%构建初始位移以及发声
u=exp(-r);% 初始振幅
t=0:0.0005:0.5;%声音最高1kHz
p=sin(2*pi*200*t);