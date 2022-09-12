%%
load('rootBesselDiff.mat')
order_max=3;%
a0=0.5; %敲击位置，(theta=0)
ca=345;a=1;L=1;
Omega=rootBesselDiff*ca/a;
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
w_nm= zeros(size(Omega));
shape_c=cell(order_max+1,order_max);% 腔的本征态
for kk=1:order_max% z方向量子数
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
end
%频率为Omega,理论上z轴需要乘sin((k+1/2)z)
%%
for nn=2:3
    for mm=1:3
        filePath=['.\pic\' num2str(nn) '.' num2str(mm)];
        mkdir(filePath)
        cnt=0;
        z_max=max(max(abs(shape_c{nn+1,mm})));
        for th=0:0.05:2*pi
            cnt=cnt+1;
            z1=shape_c{nn+1,mm}*cos(th);
            f=surfShape(z1,r,theta);
            f.EdgeColor="none";
            colormap autumn
            l=light;
            l.Position=[1,0.5,1];
            axis([-1,1,-1,1,-z_max,z_max])
            view([0.5,1,1])
            axis off
            saveas(gcf,[filePath '\' num2str(cnt) '.jpg'])
        end
    end
end