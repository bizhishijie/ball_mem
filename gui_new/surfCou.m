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
for n=1
    for m=1
        cnt=0;
        mkdir(sprintf('./couple/fig/%d.%d/',n,m));

        z_max=max(max(abs(shape_c{n+1,m})));
        omega_tmp=omega(n+1,m);
        for t=linspace(0,2*pi/omega_tmp,100)
            clf
            hold on
            cnt=cnt+1;
            s1=surfShapee(shape_c{n+1,m}*cos(omega_tmp*t)-deep,r,theta);
            s1.EdgeColor="none";
            s2=surfShapee(shape_m{n+1,m}*cos(omega_tmp*t),r,theta);
            s2.EdgeColor="none";
            l = light;
            l.Color = [1 1 1];
            l.Position = [ 0.5 -1  0.5];
            view([5,-3,6])
            colormap autumn
            xlabel('x');ylabel('y');zlabel('z');
            zlim([-z_max-deep z_max*2]);
            %             title(sprintf("%.3f",t))
            axis off
            drawnow
            %         pause(1)
            saveas(gcf,sprintf("./couple/fig/%d.%d/%d.bmp",n,m,cnt))
        end
    end
end