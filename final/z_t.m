r_num=100;
theta_num=100;% 切分的数目
order_max=10;

r=linspace(0,1,r_num)';
theta=linspace(0,2*pi-2*pi/theta_num,theta_num);
load('rootBessel.mat');
load('rootBesselDiff.mat');

load('OmegaRe.mat')


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

shape_m=cell(order_max+1,order_max);
for nn=0:order_max
    for mm=1:order_max
        z1=sqrt(2)*besselj(nn,rootBessel(nn+1,mm)*r)/...
            besselj(nn+1,rootBessel(nn+1,mm))*...
            cos(nn*theta)/sqrt(pi)/sqrt(1+(nn==0));
        % 先是r,后是theta
        shape_m{nn+1,mm}=z1;
    end
end

a=1;ca=340;
omega=Omega*ca/a;
deep=5;
clf
n_max_draw=2;
m_max_draw=2;
for n=1:n_max_draw
    for m=1:m_max_draw
        %         mkdir(sprintf('./couple/fig/%d.%d/',n,m));
        z_max=max(max(abs(shape_c{n+1,m})));
        omega_tmp=omega(n+1,m);
        z1=shape_c{n+1,m};
        z2=shape_m{n+1,m};

        z1=z1(end/2,end/2);z2=z2(end/2,end/2);
        t=linspace(0,2*pi/omega_tmp*10,100);% 默认画10个周期
        z1_t=z1*cos(omega_tmp*t);
        z2_t=z2*cos(omega_tmp*t);
        figure(1)
        subplot(n_max_draw,m_max_draw,(n-1)*m_max_draw+m)
        plot(z1_t)
        figure(2)
        subplot(n_max_draw,m_max_draw,(n-1)*m_max_draw+m)
        plot(z2_t)
    end
end