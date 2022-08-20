uiopen('.\src.fig',1)
h1=findall(gca,'type','line');
x1=get(h1,'xdata');
y1=get(h1,'ydata');
hold on

tt=40;
check=zeros(size(tt));

a =0.047;% 鼓的半径
L = 0.1;% 鼓的深度
gamma = 1.4; % 空气绝热系数
rho = 1.184;% 空气密度
% d = 3e-4; % 气球厚度
ca = 345;% 空气中声速
sigma = 0.095;% 气球面密度
pa = 1e5;% 大气压强

L_a=L/a;sigma_a_rho=sigma/a/rho;
a0=0*a;
h0=0.003;

for tii=1:length(tt)
    %%
    T=tt(tii); %张力
    cs = sqrt(T/sigma);% 橡胶中波的传播速度
    cs_ca=cs/ca;
    load('rootBessel.mat')
    load('rootBesselDiff.mat')
    for n=0:10
        fprintf("%d\t%d\n",tii,n)
        OmegaResonance(n,10, cs_ca,L_a, sigma_a_rho ,rootBessel,rootBesselDiff);
    end

    fileList=dir('./OmegaRe/*.mat');
    order_max=10;% 需要修改%%%%%%%%%%%%%%%
    Omega=zeros(length(fileList),order_max);
    for ni =1:length(fileList)
        load(strcat(fileList(ni).folder,'/', fileList(ni).name));
        while length(OmegaRe)<order_max
            OmegaRe=[OmegaRe nan];
        end
        Omega(ni,:)=OmegaRe(1:order_max);
    end
    save('OmegaRe.mat',"Omega")

    fileList=dir('./Trans/*.mat');

    % order_max=5;% 需要修改
    transMat=zeros(length(fileList),order_max,order_max);
    for ni =1:length(fileList)
        load(strcat(fileList(ni).folder,'/', fileList(ni).name));
        transMat(ni,:,:)=trans(:,1:order_max);
        % 报错就增大Omega_loop
        %     det(trans(ni));% 这个数字很接近0，计算应该是对的
    end
    save('transMat.mat',"transMat")

    fileList=dir('./Lambda/*.mat');

    % order_max=5;% 需要修改
    lambdaMat=zeros(length(fileList),order_max,order_max);
    for ni =1:length(fileList)
        load(strcat(fileList(ni).folder,'/', fileList(ni).name));
        lambdaMat(ni,:,:)=lambda;
    end
    save('lambdaMat.mat',"lambdaMat")
    %%
    load("OmegaRe.mat"); % 计算得到的OmegaRe
    load('transMat.mat'); % a^k_ni
    load("lambdaMat.mat");% lambda(n,i,j)
    load('rootBessel.mat');
    load('rootBesselDiff.mat');

    order_max=10;
    Omega=Omega(1:order_max+1,1:order_max);% 一致化
    %%
    r_num=100;
    theta_num=100;% 切分的数目
    r=linspace(0,a,r_num)';% 注意有转置，列矢*行矢比较方便 r(0)对应中间
    theta=linspace(0,2*pi,theta_num);
    % z0=zeros(r_num,theta_num);
    %%
    ab=sqrt(r.^2+a0^2-2*r*a0*cos(theta));
    cosoba=(r.^2+ab.^2-a0^2)/2./r./ab;
    coscbo=-cosoba;
    cb=(2*coscbo.*r+sqrt((2*coscbo.*r).^2-4*(r.^2-a^2)))/2;
    z0=h0*cb./(cb+ab);
    z0(isnan(z0))=h0*a/(a+a0);
    %     surfShape(z0,r,theta);
    %%
    b_mi= zeros(size(Omega));% 膜在膜的本征态上的展开系数的矩阵
    shape_m=cell(order_max+1,order_max);% 记录膜的本征态的形状的cell
    % r1=r-a/2/r_num;r2=r+a/2/r_num;
    % ds=pi*(r2.^2-r1.^2)/theta_num;
    ds=pi*2*r*a/r_num/theta_num;% 和上面两行等价
    ds=repmat(ds,1,theta_num);
    ds(:,1)=ds(:,1)/theta_num;
    for nn=0:order_max
        for mm=1:order_max
            z1=sqrt(2)*besselj(nn,rootBessel(nn+1,mm)*r/a)/...
                besselj(nn+1,rootBessel(nn+1,mm))*cos(nn*theta);
            % 先是r,后是theta  z1(end,1)约为0
            shape_m{nn+1,mm}=z1;
            b_mi(nn+1,mm)=sum(sum(z0.*z1.*ds));% 用了平方差公式，需演算
        end
    end
    % 验证计算的归一性
    shape=zeros(r_num,theta_num);
    for nn=0:order_max
        for mm=1:order_max
            shape=shape+shape_m{nn+1,mm}*b_mi(nn+1,mm);
        end
    end
    %surfShape(shape,r,theta);% 取消注释看和初态的差别
    b_mi=b_mi*h0/max(max(shape));
    %%
    b_aj=zeros(size(b_mi));
    for nn=0:order_max
        for jj=1:order_max
            b_aj(nn+1,jj)=sum(b_mi(nn+1,:)*reshape(lambdaMat(nn+1,:,jj),order_max,1));
            % 参见 OmegaResonance -13
        end
    end
    %%
    b_ck=zeros(size(b_mi));
    for nn=0:order_max
        for kk=1:order_max
            b_ck(nn+1,kk)=sum(b_aj(nn+1,:)*reshape(transMat(nn+1,kk,:),order_max,1));
            % 参见 OmegaResonance -111
        end
    end
    %%
    fn=length(x1); % 横轴划分的个数
    omega_limit=max(x1);
    x=x1;
    y=zeros(1,fn);
    for nn=0:order_max
        for mm=1:order_max
            y(ceil(Omega(nn+1,mm)*cs/a/(omega_limit/fn)))=abs(b_ck(nn+1,mm));
        end
    end
    y=y/sqrt(sum(y.^2));
%     plot(x,y(1:fn),'Color',[0 1-tii/length(tt) 1])% 需要检查y的d长度再作图
    axis([0 2000 0 1])
    drawnow
    hold on
    check(tii)=sum(y1.*y(1:fn));
end
disp(check);