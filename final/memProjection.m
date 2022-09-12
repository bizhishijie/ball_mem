load('rootBessel.mat');
load('rootBesselDiff.mat');
%%
% 定义常数 国际单位
a=0.0495-0.00216; %半径
L = 0.13-0.00216;% 鼓的深度
T=20; %张力19.9
gamma = 1.4; % 空气绝热系数
rho = 1.184;% 空气密度
ca = 346;% 空气中声速
sigma = 0.120;% 气球面密度
pa = 1e5;% 大气压强
cs = sqrt(T/sigma);% 橡胶中波的传播速度
cs_ca=cs/ca;L_a=L/a;sigma_a_rho=sigma/a/rho;
%%
order_max=20;% 需要修改
a0=0.5; %敲击位置，(theta=0)
%%
r_num=100;
theta_num=100;% 切分的数目
z_num=100;
r=linspace(0,1,r_num)';% 注意有转置，列矢*行矢比较方便 r(0)对应中间
theta=linspace(0,2*pi-2*pi/theta_num,theta_num);
z(1,1,:)=linspace(0,L,z_num);
% z0=zeros(r_num,theta_num);
% z0(floor(r_num*a0/a)+1,1)=1; % 敲击位置，方程可以后续再改
%%
% 锥形初始条件
h0=0.75;

ab=sqrt(r.^2+a0^2-2*r*a0*cos(theta));
cosoba=(r.^2+ab.^2-a0^2)/2./r./ab;
coscbo=-cosoba;
cb=(2*coscbo.*r+sqrt((2*coscbo.*r).^2-4*(r.^2-1^2)))/2;
z0=h0*cb./(cb+ab);
z0(isnan(z0))=h0*1/(1+a0);
z0=max(max(z0))-z0;
z0=exp(-z0.^2);
z0=z0-min(min(z0));
%
%%
% 球冠初始条件
% rb=0.1;% 球的半径
% zb=0.1;% 球心的高度
% z_tmp1=real(sqrt(rb^2-(r*cos(theta)-a0).^2-(r*sin(theta)).^2));
% z_tmp1(z_tmp1<0)=0;
% z0=z_tmp0+z_tmp1;
% surfShape(z0,r,theta);
%%
w_nm= zeros(order_max+1,order_max);% 膜在膜的本征态上的展开系数的矩阵
shape_m=cell(order_max+1,order_max);% 记录膜在腔的本征态下的形状的cell
% r1=r-a/2/r_num;r2=r+a/2/r_num;
% ds=pi*(r2.^2-r1.^2)/theta_num;
ds=pi*2*r/r_num/theta_num;% 和上面两行等价
ds=repmat(ds,1,theta_num);
w_cm=zeros(order_max+1,order_max);
%%
for nn=0:order_max
    for mm=1:order_max
        z1=sqrt(2)*besselj(nn,rootBessel(nn+1,mm)*r)/...
            besselj(nn+1,rootBessel(nn+1,mm))*...
            cos(nn*theta)/sqrt(pi)/sqrt(1+(nn==0));
        % 先是r,后是theta
        shape_m{nn+1,mm}=z1;
        w_cm(nn+1,mm)=sum(sum(z0.*z1.*ds));
    end
end
shape=zeros(r_num,theta_num);
for nn=0:order_max
    for mm=1:order_max
        shape=shape+shape_m{nn+1,mm}*w_cm(nn+1,mm);
    end
end
% surfShape(shape,r,theta);
% %%
Omega=rootBessel(1:(order_max+1),1:order_max);

fn=10000; % 横轴划分的个数
omega_limit=50;
x=linspace(0,omega_limit*cs/a/2/pi,fn);
y=zeros(1,fn);
for nn=0:order_max
    for mm=1:order_max
        %         y=y+(w_cm(nn+1,mm)*Omega(nn+1,mm)*cs/a/2/pi)^2./(1+(x-Omega(nn+1,mm)*cs/a/2/pi).^2/10);
        shape=shape+shape_m{nn+1,mm}*w_cm(nn+1,mm);
    end
end
% y=y/max(y);
% hold on
% plot(x,y(1:fn))% 需要检查y的长度再作图
% xlim([0 2e3])

% fn=10000; % 横轴划分的个数
% omega_limit=50;
% x=linspace(0,omega_limit*cs/a/2/pi,fn);
% y=zeros(1,fn);
% for nn=0:order_max
%     for mm=1:order_max
%         y=y+(abs(x-Omega(nn+1,mm)*cs/a/2/pi)<omega_limit*cs/a/fn/2/pi)*1;
%     end
% end
% y(y>1)=1;
% plot(x,y(1:fn))% 需要检查y的长度再作图
% xlim([0 2e3])
ii=0;
z_max=max(max(abs(shape)));
for t=0:0.00001:0.00001*200
    shape=zeros(r_num,theta_num);
    for nn=0:order_max
        for mm=1:order_max
            shape=shape+shape_m{nn+1,mm}*w_cm(nn+1,mm)*cos(Omega(nn+1,mm)*ca/a*t);
        end
    end
    ii=ii+1;
    f=surfShape(shape,r,theta);
%     colormap bone
    f.EdgeColor='none';
    axis([-1,1,-1,1,-z_max,z_max])
    l=light;
    axis off
    view(-30,70);  % 设置视点位置
    %     title(t)
    material([0.1,0.8,0.4])
    drawnow
    saveas(gcf,['./pic/' num2str(ii) '.jpg'])
    disp(t)
    pause(0.1)
end