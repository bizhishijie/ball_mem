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

            s=surfShapee(shape_m{m}*cos(omega_tmp*t),r,theta);
            s.EdgeColor="none";
            l = light;
            l.Color = [1 1 1];
            l.Position = [ 0.5 -1  0.5];
            view([5,-3,6])
            colormap autumn
            xlabel('x');ylabel('y');zlabel('z');
            zlim([-z_max(m) z_max(m)]);
            %             title(sprintf("%.3f",t))
            axis off
            drawnow
            %         pause(1)
%             saveas(gcf,sprintf("./couple/fig/%d.%d/%d.bmp",n,m,cnt))
        end
    end
end




L_a=5;
sigma_a_rho=3.3;
cs_ca=0.6;
a=1;
ca=345;cs=cs_ca*ca;

load('rootBessel.mat')
load('rootBesselDiff.mat')


parfor n=0:10
    disp(n)
    OmegaResonance(n,10, cs_ca,L_a, sigma_a_rho,rootBessel,rootBesselDiff);
end
OmegaReLoader
load('OmegaRe.mat');
rootBesselDiff=rootBesselDiff(1,:);
rootBessel=rootBessel(1,:);
Omega=Omega(1,:);
f_mem=rootBessel*ca/a/2/pi;
f_cav=rootBesselDiff*ca/a/2/pi;
f_cou=Omega*ca/a/2/pi*L_a;

f_mem=sort(f_mem(:));
f_cav=sort(f_cav(:));
f_cou=sort(f_cou(:));

clf;hold on

plot(f_cav,3*ones(size(f_cav)),'go')
plot(f_mem,2*ones(size(f_mem)),'ro')
plot(f_cou,ones(size(f_cou)),'bo')

axis([0 1000 0 4]);
set(gca,'xscale','log')