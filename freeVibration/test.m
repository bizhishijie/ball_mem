load('rootBesselDiff.mat')
order_max=5;% 需要修改!!!!这个数字太大了电脑的内存不够用
ca=345;a=1;L=1;
Omega=rootBesselDiff*ca/a;
Omega=Omega(1:order_max+1,1:order_max)/2;% 一致化
r_num=100;
theta_num=100;% 切分的数目
z_num=100;
[r,theta,z] = ndgrid(linspace(0,a,r_num),linspace(0,2*pi,theta_num),linspace(0,L,z_num));
shape_c=cell(order_max+1,order_max,order_max);% 腔的本征态
for kk=1:order_max% z方向量子数
    for nn=0:order_max
        for mm=1:order_max
            if nn==0 && mm==1
                z1=sqrt(2)*ones(size(r)).*cos(nn*theta).*...
                    sin((kk-1/2)*pi*z/L)/sqrt(pi)/sqrt(2);
            elseif nn==0
                z1=sqrt(2)*rootBesselDiff(nn+1,mm)*besselj(nn,rootBesselDiff(nn+1,mm)*r)/...
                    besselj(nn,rootBesselDiff(nn+1,mm))/sqrt(rootBesselDiff(nn+1,mm)^2-nn^2).*...
                    cos(nn*theta).*sin((kk-1/2)*pi*z/L)/sqrt(pi)/sqrt(2);
            else%mm 需要+1
                z1=sqrt(2)*rootBesselDiff(nn+1,mm+1)*besselj(nn,rootBesselDiff(nn+1,mm+1)*r)/...
                    besselj(nn,rootBesselDiff(nn+1,mm+1))/sqrt(rootBesselDiff(nn+1,mm+1)^2-nn^2).*...
                    cos(nn*theta).*sin((kk-1/2)*pi*z/L)/sqrt(pi);
            end
            % 先是r,后是theta
            shape_c{nn+1,mm,kk}=z1;
            %         w_cm(nn+1,mm)=sum(sum(z0.*z1.*ds));
        end
    end
end
tmp=shape_c{2,2,1};
contourslice(r.*cos(theta),r.*sin(theta),z,tmp,[],[],[0.5]);
