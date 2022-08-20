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
% plot(r,phi(:,3,1));% 大小下来了
sum1=sum(phi.*phi.*r,1)*r0/rDiff;% 验证归一化
% 膜的本征函数
% 
%%
%气体的压强(18)
% 第一项
% r,n,i
load('rootBesselDiff.mat','rootBesselDiff');
rDiff=1000;%将r划分的网格数
r0=1;
r=(0:r0/(rDiff-1):r0)';
n=0:10;
i=1:10;%划分网格

v=permute(rootBesselDiff,[3,1,2]);
v=repmat(v,length(r),1,1);% 重整v

sqrt2=sqrt(2);
phi1=zeros(size(v));
vr=v.*r;

for ni=n+1
    vTemp=v(:,ni,:);
    vRTemp=vr(:,ni,:);
    phi1(:,ni,:)=sqrt2*vTemp./sqrt(vTemp.^2-ni^2).*besselj(ni-1,vRTemp)./besselj(ni-1,vTemp);
end
phi1(:,1,1)=sqrt2;% (19)
clear vr v
%%
% 画关于r的函数(19)
% plot(r,phi1(:,3,1))% 数值特别大，在r=1处为0
% 对其归一化
sum1=sum(phi1.*phi1.*r,1)*r0/rDiff;% (22)
phi1=phi1./sqrt(repmat(sum1,length(r),1,1));
% plot(r,phi1(:,3,1));% 大小下来了
sum1=sum(phi1.*phi1.*r,1)*r0/rDiff;% 验证归一化
% 
%%
% 气体的压强（19）
% 第二项
% n,i,k
load('rootBesselDiff.mat','rootBesselDiff');
k=1:10;
n=0:10;
i=1:10;%划分网格

v=repmat(rootBesselDiff,1,1,length(k));% 重整v
kTemp=permute(k,[1,3,2]);
omega=sqrt((v*ca/a).^2+((repmat(kTemp,length(n),length(i),1)-1/2)*pi*ca/L));% (21)
clear kTemp
%%
% 为了(22下的式子)，将各个phi1,omega重整为4微矩阵
% r,n,i,k,      theta
phi11=repmat(phi1,1,1,1,length(k));
omega1=permute(omega,[4,1,2,3]);
omega1=repmat(omega1,size(phi11,1),1,1,1);
w=omega1.*phi11;
% 如果加上时间角度的话电脑内存不够
theta=0:0.05:2*pi;
theta=permute(theta,[5,1,4,3,2]);
whos theta
theta=repmat(theta,length(r),length(n),length(i),length(k),1);
w=repmat(w,1,1,1,1,length(theta)).*cosTheta;