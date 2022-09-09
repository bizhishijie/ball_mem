load("OmegaRe.mat"); % 计算得到的OmegaRe
load('transMat.mat'); % a^k_ni
load('rootBessel.mat');
load('rootBesselDiff.mat');
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
r_num=100;
theta_num=100;% 切分的数目
z_num=100;
r=linspace(0,1,r_num)';% 注意有转置，列矢*行矢比较方便 r(0)对应中间
theta=linspace(0,2*pi,theta_num);
%%
for nn=0:order_max
    mat_tmp=reshape(transMat(nn+1,1:order_max,1:order_max),order_max,order_max);
    for ii=1:order_max
        [~,idx]=max(abs(mat_tmp(:,ii)));
        mat_tmp(:,ii)=-mat_tmp(:,ii)*sign(mat_tmp(idx,ii));
        % 最大的数字的符号为负
    end
    [i, j] = linear_sum_assignment(mat_tmp');% 匈牙利算法
    % 行i和列j对应，由于多了一个转置，实际上是列i和行j对应
    mat_tmp=-mat_tmp(i,:);
    transMat(nn+1,1:order_max,1:order_max)=mat_tmp;
end
%%
order_max=20;
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
    end
end
shape_cou=cell(order_max,order_max);
for nn=0:order_max-1
    for mm=1:order_max
        shape_cou{nn+1,mm}=zeros(r_num,theta_num);
        for mmi=1:order_max
            shape_cou{nn+1,mm}=shape_cou{nn+1,mm}+transMat(nn+1,nn+1,mmi)*shape_c{mmi,nn+1};
        end
    end
end
m=shape_cou{2,2};
s=surfShape(m,r,theta);
s.EdgeColor="none";
% l = light;
% l.Color = [1 1 1];
view([5,-3,6])
% colormap autumn
xlabel('x');ylabel('y');zlabel('z');
zlim([min(min(m)) max(max(m))]);
axis off