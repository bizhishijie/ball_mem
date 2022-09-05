load("OmegaRe.mat"); % 计算得到的OmegaRe
load('transMat.mat'); % a^k_ni
load('rootBessel.mat');
load('rootBesselDiff.mat');
load('transMat.mat');
%%
a=0.047; %半径
L = 0.1;% 鼓的深度
gamma = 1.4; % 空气绝热系数
T=20; % 张力19.9
rho = 1.184;% 空气密度
ca = 345;% 空气中声速x
sigma = 0.120;% 气球面密度
pa = 1e5;% 大气压强
cs = sqrt(T/sigma);% 橡胶中波的传播速度
cs_ca=cs/ca;L_a=L/a;sigma_a_rho=sigma/a/rho;
%%
order_max=10;% 需要修改
a0=0; %敲击位置，(theta=0)
Omega=Omega(1:order_max+1,1:order_max)/2;% 一致化
%%
r_num=100;
theta_num=100;% 切分的数目
z_num=100;
r=linspace(0,1,r_num)';% 注意有转置，列矢*行矢比较方便 r(0)对应中间
theta=linspace(0,2*pi,theta_num);
z(1,1,:)=linspace(0,L,z_num);
% z0=zeros(r_num,theta_num);
% z0(floor(r_num*a0/a)+1,1)=1; % 敲击位置，方程可以后续再改
%%
% 锥形初始条件
h0=1;

ab=sqrt(r.^2+a0^2-2*r*a0*cos(theta));
cosoba=(r.^2+ab.^2-a0^2)/2./r./ab;
coscbo=-cosoba;
cb=(2*coscbo.*r+sqrt((2*coscbo.*r).^2-4*(r.^2-1^2)))/2;
z0=h0*cb./(cb+ab);
z0(isnan(z0))=h0*1/(1+a0);
% surfShape(z0,r,theta);
%%
% 球冠初始条件
% z0=real(sqrt(   0.5^2-(r*cos(theta)-0.5).^2-(r*sin(theta)).^2-0.2));
% z0(z0<0)=0;
% surfShape(z0,r,theta);
%%
w_nm= zeros(size(Omega));% 膜在膜的本征态上的展开系数的矩阵
% r1=r-a/2/r_num;r2=r+a/2/r_num;
% ds=pi*(r2.^2-r1.^2)/theta_num;
ds=pi*2*r/r_num/theta_num;% 和上面两行等价
ds=repmat(ds,1,theta_num);
w_cm=zeros(order_max+1,order_max);

% 计算腔的本征函数
shape_c=cell(order_max+1,order_max);
for nn=0:order_max
    for mm=1:order_max
        if nn==0 && mm==1
            z1=sqrt(2)*ones(size(r))*cos(nn*theta)/sqrt(pi)/sqrt(2);
        elseif nn==0
            z1=sqrt(2)*rootBesselDiff(nn+1,mm)*besselj(nn,rootBesselDiff(nn+1,mm)*r)/...
                besselj(nn,rootBesselDiff(nn+1,mm))/sqrt(rootBesselDiff(nn+1,mm)^2-nn^2)*...
                cos(nn*theta)/sqrt(pi)/sqrt(2);
        else%mm 需要+1
            z1=sqrt(2)*rootBesselDiff(nn+1,mm+1)*besselj(nn,rootBesselDiff(nn+1,mm+1)*r)/...
                besselj(nn,rootBesselDiff(nn+1,mm+1))/sqrt(rootBesselDiff(nn+1,mm+1)^2-nn^2)*...
                cos(nn*theta)/sqrt(pi);
        end
        % 先是r,后是theta
        shape_c{nn+1,mm}=z1;
        w_cm(nn+1,mm)=sum(sum(z0.*z1.*ds));
    end
end
% 计算膜的本征函数
shape_m=cell(order_max,order_max);
for nn=0:order_max
    for mm=1:order_max
        z1=sqrt(2)*besselj(nn,rootBessel(nn+1,mm)*r)/...
            besselj(nn+1,rootBessel(nn+1,mm))*...
            cos(nn*theta)/sqrt(pi)/sqrt(1+(nn==0));
        % 先是r,后是theta
        shape_m{nn+1,mm}=z1;
    end
end
shape=zeros(r_num,theta_num);
for nn=0:order_max
    for mm=1:order_max
        shape=shape+shape_c{nn+1,mm}*w_cm(nn+1,mm);
    end
end
%膜的初态按照腔的本征态展开
% surfShape(shape,r,theta);% 取消注释看和初态的差别
%%
% 角标应为n m k
w_cn=zeros(size(w_cm));
for nn=0:order_max
    %     w_cn(nn+1,:)=((reshape(transMat(nn+1,:,:),order_max,order_max))^-1*w_cm(nn+1,1:order_max)')';
    w_cn(nn+1,:)=w_cm(nn+1,1:order_max)*reshape(transMat(nn+1,1:order_max,1:order_max),order_max,order_max);
end
%%
% 计算耦合体系
z_1=[];
figure(1)
z0=0.5;

for t=0:0.0005:0.1
    shape=zeros(r_num,theta_num);
    for nn=0:order_max
        for mm=1:order_max
            shape=shape+shape_c{nn+1,mm}*w_cn(nn+1,mm)*cos(Omega(nn+1,mm)*ca/a*t)*sin((k-1/2)*pi*z0);
        end
    end
    surfShape(shape,r,theta);
    pause(0.01)
    z_tmp=shape(end/2,end/2);
    z_1=[z_1 z_tmp];
end

figure(2)
plot(z_1)
title('耦合体系腔的振动')
z_1=[];

figure(3)
for t=0:0.0005:0.1
    shape=zeros(r_num,theta_num);
    for nn=0:order_max
        for mm=1:order_max
            shape=shape+shape_m{nn+1,mm}*w_cn(nn+1,mm)*cos(Omega(nn+1,mm)*ca/a*t);
        end
    end
    surfShape(shape,r,theta)
    pause(0.01)
    z_tmp=shape(end/2,end/2);
    z_1=[z_1 z_tmp];
end

figure(4)
plot(z_1)
title('耦合体系膜的振动')
figure(5)
surfShape(shape_c{1,1},r,theta);% 后面的是两个量子数，n+1,m
figure(6)
surfShape(shape_m{1,1},r,theta)% 后面的是两个量子数，n+1,m