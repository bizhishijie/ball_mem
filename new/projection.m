%%
% 定义常数 国际单位
a = 0.2;% 鼓的半径
ca = 340;% 空气中声速
cs = 1;% 橡胶中波的传播速度
L = 0.3;% 鼓的深度
gamma = 1.4; % 空气绝热系数
rho = 1.5e3;%橡胶密度
d = 1e-4; % 气球厚度
sigma = rho*d;% 气球面密度
pa = 1e5;% 大气压强
r0=0.1;% 敲击位置
% 修改参数需要重新运行 OmegaReSaver
load('rootBessel.mat')
load('rootBesselDiff.mat')
load('OmegaRe.mat')
% 以上三个文件可能部分需要重写
orderN=12;
orderMax=30;% ijm的阶数
%需要自行检验和OmegaRe中的阶数是否相等
% 获得鼓的本征频率
%%
% 计算正交基
syms r z theta
% syms r_a z_L
r_a=r/a;
z_L=z/L;%归一化
rootBessel=rootBessel(1:orderN+1,1:orderMax);
rootBesselDiff=rootBesselDiff(1:orderN+1,1:orderMax);
Omega=Omega(1:orderN+1,1:orderMax)*ca/a;
n = (0:orderN)';m = 1:orderMax;
n = repmat(n,1,orderMax);
m = repmat(m, orderN);
phi_n_m=sqrt(2)*besselj(n, rootBessel*r_a)./besselj(n+1, rootBessel);
% phi部分完毕

epsilon=L/a*sqrt(abs((Omega*a/ca).^2-rootBesselDiff.^2));
zeta=(Omega*a/ca<rootBesselDiff).*(cosh(epsilon-epsilon*z_L)./epsilon)+...
    (Omega*a/ca==rootBesselDiff)+...
    (Omega*a/ca>rootBesselDiff).*(cos(epsilon*z_L)+tan(epsilon).*sin(epsilon*z_L));
% zeta部分完毕
w=phi_n_m.*zeta.*cos(n*theta);
% 得到正交基
% 将其归一，r在=0~1，z在=0~1
heft=eval(int(int(int(w.^2*r*pa,r,[0,a]),z,[0,L]),theta ,[0,2*pi]));
w=w./heft;

%%
% 根据初态计算分量
d1=sqrt(r0^2+r^2+2*r*r0*cos(theta));
d2=a-r;
p_n_m_k=int(int(int(w*d1/(d1+d2),r,[0,a]),z,[0,L]),theta ,[0,2*pi]) ;
pa_n_m_k=eval(p_n_m_k);
%%
syms t
syms r z theta

n = (0:orderN)';m = 1:orderMax;
n = repmat(n,1,orderMax);
m = repmat(m, orderN);
phi_n_m=sqrt(2)*besselj(n, rootBessel*r_a)./besselj(n+1, rootBessel);
% phi部分完毕
epsilon=L/a*sqrt(abs((Omega*a/ca).^2-rootBesselDiff.^2));
zeta=(Omega*a/ca<rootBesselDiff).*(cosh(epsilon-epsilon*z_L)./epsilon)+...
    (Omega*a/ca==rootBesselDiff)+...
    (Omega*a/ca>rootBesselDiff).*(cos(epsilon*z_L)+tan(epsilon).*sin(epsilon*z_L));
% zeta部分完毕
p=pa_n_m_k.*phi_n_m.*zeta.*cos(n*theta).*sin(Omega*t);
w_t2=-diff(p,z)/rho/L;
w_t2=subs(w_t2,z,0);
w=int(int(w_t2,t),t);
w_sum=sum(sum(w));
% save('main.mat')
%%
% 可视化
rSize=100;thetaSize=30;
r_loop=linspace(0,a,rSize);
theta_loop=linspace(0,2*pi,thetaSize);
r_loop_1=repmat(r_loop,thetaSize,1);
theta_loop_1=repmat(theta_loop',1,rSize);

z=zeros(rSize,thetaSize);
t0=1;t_diff=1e-4;t1=1.1;
for time=t0:t_diff:t1
    w_temp=subs(w_sum,t,time);
    for ii=1:rSize
        temp=subs(w_temp,r,r_loop(ii));
        for jj=1:thetaSize
            %         z(ii,jj)=subs(temp,theta,theta_loop(jj));
            z(ii,jj)=temp;
        end
    end

    surf(r_loop_1.*cos(theta_loop_1),r_loop_1.*sin(theta_loop_1),z')
    title(time);
    %     axis([-a,a,-a,a,-5e-7,5e-5]);
    view(-30,70);  % 设置视点位置
    drawnow

    saveas(gcf,['./pic1/' num2str((time-t0)/t_diff) '.jpg'])
    disp(time)
end