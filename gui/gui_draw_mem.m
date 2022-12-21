function gui_draw_mem
th_n= app.EditField_theta.Value;
r_n=app.EditField_r.Value;
load('.\para.mat')
load('.\rootBessel.mat');
load('.\rootBesselDiff.mat');
load('.\OmegaRe.mat')
%a=1; %归一化后
r_num=100;
theta_num=100;% 切分的数目
order_max=10;
r=linspace(0,1,r_num)';
theta=linspace(0,2*pi-2*pi/theta_num,theta_num);
n0=th_n;

shape_m=cell(1,order_max);
for mm=1:order_max
    z1=sqrt(2)*besselj(n0,rootBessel(n0+1,mm)*r)/...
        besselj(n0+1,rootBessel(n0+1,mm))*...
        cos(n0*theta)/sqrt(pi)/sqrt(1+(n0==0));
    % 先是r,后是theta
    shape_m{mm}=z1;
end

z_max=zeros(1,order_max);
for mm=1:order_max
    z_max(mm)=max(max(shape_m{mm}));
end

n=th_n;
m=r_n;
cnt=0;
mkdir(sprintf('./couple/fig/%d.%d/',n,m));
omega_tmp=Omega(n+1,m);

for t=linspace(0,2*pi/omega_tmp,100)
    cnt=cnt+1;
    s=surfShape(app.UIAxes2_2,shape_m{m}*cos(omega_tmp*t),r,theta);
    s.EdgeColor="none";
    axis(app.UIAxes2_2,[-1,1,-1,1,-1.5,1.5]);
    load('s_w0.mat');
    load('s_h0.mat');
    view(app.UIAxes2_2,see_where0,see_high0);  % 设置视点位置
    axis off
    drawnow
    pause(0.05)
    %        saveas(gcf,sprintf("./couple/fig/%d.%d/%d.bmp",n,m,cnt))
end
end