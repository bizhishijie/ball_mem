load("OmegaRe.mat"); % 计算得到的OmegaRe
load('transMat.mat'); % a^k_ni
load('rootBessel.mat');
load('rootBesselDiff.mat');
%%
a=0.047; %半径
L = 0.1;% 鼓的深度
gamma = 1.4; % 空气绝热系数
T=19; % 张力19.9
rho = 1.184;% 空气密度
ca = 345;% 空气中声速
sigma = 0.120;% 气球面密度
pa = 1e5;% 大气压强
cs = sqrt(T/sigma);% 橡胶中波的传播速度
cs_ca=cs/ca;L_a=L/a;sigma_a_rho=sigma/a/rho;
%%
order_max=20;% 需要修改
a0=1; %敲击位置，(theta=0)
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
z0=max(max(z0))-z0;
z0=exp(-z0.^2);
z0=z0-min(min(z0));
% surfShape(z0,r,theta);
%%
% 球冠初始条件
% z0=real(sqrt(0.5^2-(r*cos(theta)-0.5).^2-(r*sin(theta)).^2-0.2));
% z0(z0<0)=0;
% surfShape(z0,r,theta);
%%
w_nm= zeros(size(Omega));% 膜在膜的本征态上的展开系数的矩阵
shape_c=cell(order_max+1,order_max);% 记录膜在腔的本征态下的形状的cell
% r1=r-a/2/r_num;r2=r+a/2/r_num;
% ds=pi*(r2.^2-r1.^2)/theta_num;
ds=pi*2*(r+a/r_num)/r_num/theta_num;% 和上面两行等价
ds=repmat(ds,1,theta_num);
w_cm=zeros(order_max+1,order_max);
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
shape=zeros(r_num,theta_num);
for nn=0:order_max
    for mm=1:order_max
        shape=shape+shape_c{nn+1,mm}*w_cm(nn+1,mm);
    end
end
%膜的初态按照腔的本征态展开
% surfShape(shape,r,theta);% 取消注释看和初态的差别
%%
% for nn=0:order_max
%     mat_tmp=reshape(transMat(nn+1,1:order_max,1:order_max),order_max,order_max);
%     for ii=1:order_max
%         [~,idx]=max(abs(mat_tmp(:,ii)));
%         mat_tmp(:,ii)=-mat_tmp(:,ii)*sign(mat_tmp(idx,ii));
%     end
%     [i, j] = linear_sum_assignment(mat_tmp');
%     % 行i和列j对应，由于多了一个转置，实际上是列i和行j对应
%     mat_tmp=-mat_tmp;
%     transMat(nn+1,1:order_max,1:order_max)=mat_tmp;
% end
%%
% 角标应为n m k
w_cn=zeros(size(w_cm));
for nn=0:order_max
    %     w_cn(nn+1,:)=((reshape(transMat(nn+1,:,:),order_max,order_max))^-1*w_cm(nn+1,1:order_max)')';
    w_cn(nn+1,:)=w_cm(nn+1,1:order_max)*reshape(transMat(nn+1,1:order_max,1:order_max),order_max,order_max);
end
%%
% fn=10000; % 横轴划分的个数
% omega_limit=10;
% x=linspace(0,omega_limit*ca/a/2/pi,fn);
% y=zeros(1,fn);
% for nn=0:order_max
%     for mm=1:order_max
%         y=y+abs(w_cn(nn+1,mm)*Omega(nn+1,mm))./(1+(x-Omega(nn+1,mm)*ca/a/2/pi).^2/10);
%     end
% end
% y=y/max(y);
% hold on
% plot(x,y(1:fn))% 需要检查y的长度再作图
% xlim([0 2e3])

% fn=10000; % 横轴划分的个数
% omega_limit=20;
% x=linspace(0,omega_limit*ca/a/2/pi,fn);
% y=zeros(1,fn);
% for nn=0:order_max
%     for mm=1:order_max
%         y=y+(abs(x-Omega(nn+1,mm)*ca/a/2/pi)<omega_limit*ca/a/fn/2/pi)*1;
%     end
% end
% y(y>1)=1;
% plot(x,y(1:fn))% 需要检查y的长度再作图
% xlim([0 800])
%%
ii=1;
mkdir('./pic')
for t=0:0.0001:0.02
    shape=zeros(r_num,theta_num);
    for nn=0:order_max
        for mm=1:order_max
            shape=shape+shape_c{nn+1,mm}*w_cn(nn+1,mm)*cos(Omega(nn+1,mm)*ca/a*t);
        end
    end
    surfShape(shape,r,theta);
    view(-30,70);  % 设置视点位置
    %     axis([-a,a,-a,a,-0.003,0]); % 显示不全调节最后两个参数
    title(t)
    drawnow
    %     saveas(gcf,['./pic/' num2str(ii) '.jpg'])
    disp(t)
    ii=ii+1;
    pause(0.1)
end