load("OmegaRe.mat");
cs_ca=0.3;sigma_a_rho=3.3;L_a=1.5; % 论文中B鼓
order_max=10;

% 定义常数 国际单位
a = 1;% 鼓的半径
ca = 346;% 空气中声速
cs = cs_ca*ca;% 橡胶中波的传播速度
L = L_a*a;% 鼓的深度
gamma = 1.4; % 空气绝热系数
rho = 1.184;% 空气密度
% d = 3e-4; % 气球厚度
sigma = sigma_a_rho*a*rho;% 气球面密度
pa = 1e5;% 大气压强
a0=0.2*a; %敲击位置，(theta=0)

Omega=Omega(1:order_max+1,1:order_max);% 一致化
%%
load('rootBessel.mat')
load('rootBesselDiff.mat')

r_num=100;
theta_num=100;% 切分的数目
r=linspace(0,a,r_num)';% 注意有转置，列矢*行矢比较方便 r(0)对应中间
theta=linspace(0,2*pi,theta_num);
z0=zeros(r_num,theta_num);
z0(floor(r_num*a0/a)+1,1)=1; % 敲击位置，方程可以后续再改

factor_mem= zeros(size(Omega));% 膜在膜的本征态上的展开系数的矩阵
shape_mem=cell(order_max+1,order_max);% 记录本征态的形状的cell
for nn=0:order_max
    for mm=1:order_max
        if ~isnan(Omega(nn+1,mm))
            z1=sqrt(2)*besselj(nn,rootBessel(nn+1,mm)*r/a)/besselj(nn+1,rootBessel(nn+1,mm))*...
                cos(nn*theta);% 先是r,后是theta  z1(end,1)约为0
            shape_mem{nn+1,mm}=z1;
            %             factor_mem(nn+1,mm)=sum(sum(z0.*z1))*a/r_num*2*pi/theta_num;% 后面的系数可能不对
            %         之后这里可以修改后面的微元的面积使得其符合数学
            factor_mem(nn+1,mm)=sum(sum(z0.*z1*2.*r*a/r_num*theta_num));% 用了平方差公式，需演算
        end
    end
end
factor_mem(factor_mem<1e-5)=0;% 振幅太小的不必计算，这行尽量不要删去，数字太小的话计算有误差

% factor_mem(nn,mm)表示膜在膜的本征态上的展开系数
%%
load("trans.mat"); % load变换矩阵
factor_cav=zeros(size(factor_mem));
for ni=1:order_max+1
    factor_cav(ni,:)=factor_mem(ni,:)*trans{ni};
end% 获得了在腔中的展开系数

load('Ktrans.mat');
factor_mem_new=zeros(size(factor_mem));
% Ktrans(order,:,:)表示第ii阶的K矩阵
for ni=0:order_max
    K=Ktrans{ni+1};
    for mi=1:order_max
        K_tmp=reshape(K(mi,:,:),order_max,order_max);
        factor_mem_new=factor_mem_new+(K_tmp^-1*factor_cav')';
        % 数字算的很大，达到了几十的量级，怀疑
    end
end
factor_mem_new=factor_mem_new/order_max^2;
%%
% 画频谱
% figure
fn=1000; % 横轴划分的个数
x=linspace(0,10*cs/a,fn);
y=zeros(1,fn);
for nn=0:order_max
    for mm=1:order_max
        if ~isnan(factor_cav(nn+1,mm))% 我不清楚应该是哪个系数
            y(floor(Omega(nn+1,mm)/(10/fn)))=abs(factor_cav(nn+1,mm));
        end
    end
end
plot(x,y(1:fn))% 需要检查y的长度再作图
%%
% 画空间形状，取消下面注释可以保存，之后使用createGIF可以画动图
for t=0:0.001:0.1
    shape=zeros(r_num,theta_num);
    for nn=0:order_max
        for mm=1:order_max
            if ~isnan(factor_mem_new(nn+1,mm))
                shape=shape+shape_mem{nn+1,mm}*factor_mem_new(nn+1,mm)*cos(rootBessel(nn+1,mm)*cs/a*t);
            end
        end
    end
    surf(r*cos(theta),r*sin(theta),shape)
    view(-30,70);  % 设置视点位置
    axis([-1,1,-1,1,-70,70])
    drawnow
    %     saveas(gcf,['./pic/' num2str(t) '.jpg'])
    disp(t)
    pause(0.1)
end
