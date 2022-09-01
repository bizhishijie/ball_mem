r_num=100;
theta_num=100;% 切分的数目
order_max=10;

r=linspace(0,1,r_num)';
theta=linspace(0,2*pi-2*pi/theta_num,theta_num);
load('rootBessel.mat');
load('rootBesselDiff.mat');
load('OmegaRe.mat')

n0=1;
nn=n0;

shape_m=cell(1,order_max);
for mm=1:order_max
    z1=sqrt(2)*besselj(nn,rootBessel(nn+1,mm)*r)/...
        besselj(nn+1,rootBessel(nn+1,mm))*...
        cos(nn*theta)/sqrt(pi)/sqrt(1+(nn==0));
    % 先是r,后是theta
    shape_m{mm}=z1;
end

a=1;ca=340;
omega=Omega*ca/a;
deep=5;
z_max=zeros(1,order_max);
for mm=1:order_max
    z_max(mm)=max(max(shape_m{mm}));
end
for n=n0
    for m=1:4
        cnt=0;
        mkdir(sprintf('./couple/fig/%d.%d/',n,m));
        omega_tmp=omega(n+1,m);     
        for t=linspace(0,2*pi/omega_tmp,100)
            clf
            hold on
            cnt=cnt+1;

            s=surfShape(shape_m{m}*cos(omega_tmp*t),r,theta);
            s.EdgeColor="none";
            l = light;
            l.Color = [1 1 1];z
            view([5,-3,6])
            colormap autumn
            xlabel('x');ylabel('y');zlabel('z');
            zlim([-z_max(m) z_max(m)]);
            %             title(sprintf("%.3f",t))
            axis off
            drawnow
            %         pause(1)
            saveas(gcf,sprintf("./couple/fig/%d.%d/%d.bmp",n,m,cnt))
        end
    end
end