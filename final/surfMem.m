r_num=100;
theta_num=100;% 切分的数目
order_max=10;

r=linspace(0,1,r_num)';
theta=linspace(0,2*pi-2*pi/theta_num,theta_num);
load('rootBessel.mat');
load('rootBesselDiff.mat');


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
a=1;cs=20;
omega=rootBessel*cs/a;
for n=0:2
    for m=1:3
        cnt=0;
        mkdir(sprintf('./mem/%d.%d/',n,m));
        omega_tmp=omega(n+1,m);
        z_max=max(max(abs(shape_m{n+1,m})));
        for t=linspace(0,2*pi/omega_tmp,100)
            cnt=cnt+1;
            s=surfShape(shape_m{n+1,m}*cos(omega_tmp*t),r,theta);
            s.EdgeColor="none";
            l = light;
            l.Color = [1 1 1];
            l.Position = [ 0.5 -1  0.5];
            view([5,-3,6])
            colormap autumn
            xlabel('x');ylabel('y');zlabel('z');
            zlim([-z_max z_max]);
            %             title(sprintf("%.3f",t))
            axis off
            drawnow
            %         pause(1)
            saveas(gcf,sprintf("./mem/%d.%d/%d.bmp",n,m,cnt))
        end
    end
end