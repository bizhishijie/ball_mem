%%
load('rootBessel.mat')
order_max=10;% 需要修改
a0=0.5; %敲击位置，(theta=0)
cs=20;a=1;L=1;
Omega=rootBessel*cs/a;
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
%%
w_nm= zeros(size(Omega));% 膜在膜的本征态上的展开系数的矩阵
shape_m=cell(order_max+1,order_max);% 记录膜在腔的本征态下的形状的cell
% r1=r-a/2/r_num;r2=r+a/2/r_num;
% ds=pi*(r2.^2-r1.^2)/theta_num;
ds=pi*2*r/r_num/theta_num;% 和上面两行等价
ds=repmat(ds,1,theta_num);
w_mm=zeros(order_max+1,order_max);
for nn=0:order_max
    for mm=1:order_max
        z1=sqrt(2)*besselj(nn,rootBessel(nn+1,mm)*r)/...
            besselj(nn+1,rootBessel(nn+1,mm))*...
            cos(nn*theta)/sqrt(pi)/sqrt(2);
        % 先是r,后是theta
        shape_m{nn+1,mm}=z1;
        w_mm(nn+1,mm)=sum(sum(z0.*z1.*ds));% 注释本行
    end
end
%频率为Omega
%%
ii=1;
mkdir('./pic')
for t=0:0.0001:0.02
    shape=zeros(r_num,theta_num);
    for nn=0:order_max
        for mm=1:order_max
            shape=shape+shape_m{nn+1,mm}*w_mm(nn+1,mm)*cos(Omega(nn+1,mm)*cs/a*t);
        end
    end
    surfShape(shape,r,theta)
    view(-30,70);  % 设置视点位置
    %     axis([-a,a,-a,a,-0.003,0]); % 显示不全调节最后两个参数
    title(t)
    drawnow
    %     saveas(gcf,['./pic/' num2str(ii) '.jpg'])
    disp(t)
    ii=ii+1;
    pause(0.1)
end