%%
% 定义常数 国际单位
a = 0.2;% 鼓的半径
ca = 340;% 空气中声速
cs = 1800;% 橡胶中声速
L = 0.3;% 鼓的深度
gamma = 1.4; % 空气绝热系数
ro = 1.5e3;%橡胶密度
d = 1e-4; % 气球厚度
sigma = ro*d;% 气球面密度
pa = 1e5;% 大气压强
%%
syms t r theta
% 求解时间 0~1 0~2pi
load('rootBessel.mat');
load('rootBesselDiff.mat');
% 加载贝塞尔及其极值点
%%
% 计算二维圆膜振动的方程的解
% n, m
% syms w
nOrder = 50;mOrder = 50;
n = 0:nOrder;m = 1:mOrder;
% 贝塞尔根的阶数
n = repmat(n, mOrder, 1)';
m = repmat(m, nOrder+1, 1);
% 构建本征函数集
w = sqrt2*besselj(n, rootBessel*r)./besselj(n+1, rootBessel).*cos(n*theta).*sin(rootBessel*cs/a*t);
omeganm = rootBessel*cs/a;
% 包含归一化
%%
% 膜受外界驱动力受迫振动
% syms pb pe omega
syms omega
% omega = 2*pi*200;% 认为驱动频率是200Hz
psnm = pa;
qsnm = pa;
% 定义压强的量级，认为是一个大气压
pb = psnm*besselj(n, rootBessel*r)./besselj(n+1, rootBessel).*cos(n*theta).*sin(omega*t);
qb = qsnm*besselj(n, rootBessel*r)./besselj(n+1, rootBessel).*cos(n*theta).*sin(omega*t);
wsnm = 1/sigma./(omega^2-omeganm.^2)*(psnm-qsnm);
% 得到受迫振动振幅
% 缺少归一化
%%
% 计算空气腔自身振动
% n, i, k
% syms p z z为归一化后的纵向坐标
syms z
nOrder = 50;iOrder = 50;
n = 0:nOrder;i = 1:iOrder;
kOrder = 60;
k = 1:kOrder;
% 贝塞尔根的阶数 以及纵向驻波的阶数
n = repmat(n, iOrder, 1, kOrder);n = permute(n, [2, 1, 3]);
i = repmat(i, kOrder, 1, nOrder+1);i = permute(i, [3, 2, 1]);
k = repmat(k, nOrder+1, 1, iOrder);k = permute(k, [1, 3, 2]);
% 按理不需要后面的重整维度，需优化
% 构建本征函数集
rootBesselDiffTemp = repmat(rootBesselDiff, 1, 1, kOrder);
omegaanmk = ca*sqrt((rootBesselDiffTemp/a).^2+((k-1/2)*pi/L).^2);
p = sqrt2*besselj(n, rootBesselDiffTemp*r)./sqrt(rootBesselDiffTemp.^2-n.^2)./besselj(n, rootBesselDiffTemp).*sin((k-0.5)*pi*z).*cos(n*theta).*sin(omegaanmk*t);
p(1, 1) = sqrt2;
save('main.mat');
%%
% 空气腔受迫振动
syms zeta wni
Omega = a*omega/ca;% 大写OMEGA

nOrder = 50;mOrder = 50;
n = 0:nOrder;i = 1:mOrder;
n = repmat(n, mOrder, 1)';
m = repmat(m, nOrder+1, 1);

% omega需要较大才能体现出区别
epsl = L*(sqrt(rootBesselDiff.^2-Omega^2)/a);
zeta = (Omega<rootBesselDiff).*cosh(epsl-epsl*z)./cosh(epsl)+...
    (Omega == rootBesselDiff)*1+...
    (Omega>rootBesselDiff).*(cos(epsl*z)+tan(epsl).*sin(epsl));
p = zeta.*besselj(n, rootBesselDiff*r).*cos(n*theta).*sin(omega*t);
%%
% syms Omega
% n, i, j
% n, i,m
nOrder = 50;iOrder = 50;jOrder = 50;mOrder=50;
n = 0:nOrder;i = 1:iOrder;j = 1:jOrder;m=1:mOrder;

% 贝塞尔根的阶数 以及纵向驻波的阶数
n = repmat(n, iOrder, 1, jOrder);n = permute(n, [2, 1, 3]);
i = repmat(i, jOrder, 1, nOrder+1);i = permute(i, [3, 2, 1]);
j = repmat(j, nOrder+1, 1, iOrder);j = permute(j, [1, 3, 2]);
rootBesselDiffTemp = repmat(rootBesselDiff, 1, 1, jOrder);
rootBesselTemp=repmat(rootBessel, 1, 1, jOrder);

for Omega=0:0.01:2*pi*500

    A=zeros(nOrder,iOrder,jOrder);
    lambda=zeros(nOrder,iOrder,mOrder);
    lambda=rootBesselDiffTemp./sqrt(rootBesselDiffTemp.^2-n.^2)*2.*rootBesselTemp./(rootBesselDiffTemp.^2-rootBesselTemp.^2);
    lambda(1,:,:)=2*rootBesselTemp(1,:,:)./(rootBesselDiffTemp(1,:,:).^2-rootBesselTemp(1,:,:).^2);

    lambda1=lambda;% n,i,m
    lambda2=lambda;% n,j,m
    % A=a^2*ro/gamma/pa/sigma.*sum(lambda./((cs/ca)^2*rootBesselTemp(:,ii,ii).^2)-Omega^2,3);
    % 顶不住了，需要for了
    for ni=(0:nOrder)+1
        for ii=1:iOrder
            for ji=1:jOrder
                sumTemp=0;
                for mi=1:mOrder
                    sumTemp=sumTemp+lambda1(ni,ii,mi)*lambda(ni,ji,mi)/((cs/ca)^2*rootBessel(ni,mi)-Omega^2);
                end
                A(ni,ii,ji)=a^2*ro/gamma/pa/sigma*sumTemp;
            end
        end
    end

    K=zeros(nOrder,iOrder,jOrder);
    % K (:,i==j,j==i) = (Omega<rootBesselDiffTemp)*gamma*pa/a*Omega^2./sqrt(rootBesselDiffTemp.^2-Omega^2).*coth(L/a.*sqrt(rootBesselDiffTemp.^2-Omega^2))+...
    %     (Omega>rootBesselDiffTemp)*gamma*pa/a*Omega^2./sqrt(rootBesselDiffTemp.^2-Omega^2).*coth(L/a.*sqrt(rootBesselDiffTemp.^2-Omega^2));
    % K (rootBesselDiffTemp == Omega,i==j,j==i)=inf;
    % 既然都用了一个for for for for 了，为啥不摆烂捏

    for ni=(0:nOrder)+1
        for ii=1:iOrder
            for ji=1:jOrder
                if ii~=ji
                    continue
                elseif Omega<rootBessel(ni,ii)
                    K(ni,ii,ji)=-gamma*pa/a*Omega^2/sqrt(rootBesselDiff(ni,ii)^2-Omega^2)*coth(L/a*sqrt(rootBesselDiff(ni,ii)^2-Omega^2));
                elseif Omega==rootBessel(ni,ii)
                    K(ni,ii,ji)=inf;
                elseif Omega>rootBessel(ni,ii)
                    K(ni,ii,ji)=gamma*pa/a*Omega^2/sqrt(Omega^2-rootBesselDiff(ni,ii)^2)*coth(L/a*sqrt(Omega^2-rootBesselDiff(ni,ii)^2));
                end
            end
        end
    end
    mat=1+K.*A;
    for ni =(0:nOrder)+1
        tempMat=reshape(mat(ni,:,:),iOrder,jOrder);
        tempDet=det(tempMat);
        if tempDet<1e-10
            fprintf("%d\t%.2f\t%f\n",ni,Omega,tempDet);
        end
    end
end