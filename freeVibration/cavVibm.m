%%
load('rootBesselDiff.mat')
order_max=5;% 需要修改!!!!这个数字太大了电脑的内存不够用
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
% 锥形初始条件
h0=1;

ab=sqrt(r.^2+a0^2-2*r*a0*cos(theta));
cosoba=(r.^2+ab.^2-a0^2)/2./r./ab;
coscbo=-cosoba;
cb=(2*coscbo.*r+sqrt((2*coscbo.*r).^2-4*(r.^2-1^2)))/2;
z0=h0*cb./(cb+ab);
z0(isnan(z0))=h0*1/(1+a0);
z0=repmat(z0,1,1,z_num);
%%
w_nm= zeros(size(Omega));
shape_c=cell(order_max+1,order_max,order_max);% 腔的本征态
% r1=r-a/2/r_num;r2=r+a/2/r_num;
% ds=pi*(r2.^2-r1.^2)/theta_num;
dv=pi*2*r/r_num/theta_num*L/z_num;% 和上面两行等价
dv=repmat(dv,1,theta_num);
% w_cm=zeros(order_max+1,order_max);
for kk=1:order_max% z方向量子数
    for nn=0:order_max
        for mm=1:order_max
            if nn==0 && mm==1
                z1=repmat(sqrt(2)*ones(size(r))*cos(nn*theta),1,1,z_num).*...
                    repmat(sin((kk-1/2)*pi*z/L),r_num,theta_num,1)/sqrt(pi)/sqrt(2);
            elseif nn==0
                z1=repmat(sqrt(2)*rootBesselDiff(nn+1,mm)*besselj(nn,rootBesselDiff(nn+1,mm)*r)/...
                    besselj(nn,rootBesselDiff(nn+1,mm))/sqrt(rootBesselDiff(nn+1,mm)^2-nn^2)*...
                    cos(nn*theta),1,1,z_num).*...
                    repmat(sin((kk-1/2)*pi*z/L),r_num,theta_num,1)/sqrt(pi)/sqrt(2);
            else%mm 需要+1
                z1=repmat(sqrt(2)*rootBesselDiff(nn+1,mm+1)*besselj(nn,rootBesselDiff(nn+1,mm+1)*r)/...
                    besselj(nn,rootBesselDiff(nn+1,mm+1))/sqrt(rootBesselDiff(nn+1,mm+1)^2-nn^2)*...
                    cos(nn*theta),1,1,z_num).*...
                    repmat(sin((kk-1/2)*pi*z/L),r_num,theta_num,1)/sqrt(pi);
            end
            % 先是r,后是theta
            shape_c{nn+1,mm,kk}=z1;
            w_cm(nn+1,mm)=sum(sum(sum(z0.*z1.*dv)));
        end
    end
end
%频率为Omega,理论上z轴需要乘sin((k+1/2)z)
%%
tmp=shape_c{2,2,1};
figure(1)
contourfShape(tmp,r,theta,z,0.5);
