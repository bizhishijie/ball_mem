load('rootBessel.mat');
load('rootBesselDiff.mat');
order_max=3;
shape_m=cell(order_max,order_max);

r_num=100;
theta_num=100;% 切分的数目
z_num=100;
r=linspace(0,1,r_num)';% 注意有转置，列矢*行矢比较方便 r(0)对应中间
theta=linspace(0,2*pi,theta_num);

nn=1;mm=2;
z1=sqrt(2)*besselj(nn,rootBessel(nn+1,mm)*r)/...
    besselj(nn+1,rootBessel(nn+1,mm))*...
    cos(nn*theta)/sqrt(pi)/sqrt(1+(nn==0));
% 先是r,后是theta
shape_m{nn+1,mm}=z1;
%         subplot(order_max,order_max,nn*order_max+mm)
f=surfShape(shape_m{nn+1,mm},r,theta);
f.EdgeColor='none';
axis([-1,1,-1,1,-5,5])
title(['[' num2str(nn) ',' num2str(mm) ']'],'Position',[0.7,1],'FontSize',15)
axis off
set(gca,'CameraViewAngle',8)

set(gca,'color',[1 1 1])
%%
% clf
% hold on
% x=linspace(0,10,500);
% cString = {'#93B24D', '#E6E294', '#C37F95', '#79C4D1', '#38375A', ...
%     '#48948B', '#135A67', '#69356E', '#9C3F2F', '#C4D9A7'};
%
% for nn=0:5
%     y=besselj(nn,x);
%     plot(x,y,'Color',cString{nn+1},'LineWidth',3);
% end
% legend('J_0','J_1','J_2','J_3','J_4','J_5')
% set(gca,'color',[1 1 1])
% ylabel('J_n(x)','FontSize',15)
% xlabel('x','FontSize',15)
% set(gca,'FontSize',12   );
% grid on
%%
shape_c=cell(order_max,order_max);
% nn=1;mm=2;
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
        shape_c{nn+1,mm}=z1;
    end
end
% 先是r,后是theta
shape_c{nn+1,mm}=z1;
% subplot(order_max,order_max,nn*order_max+mm)
f=surfShape(shape_c{nn+1,mm},r,theta);
f.EdgeColor='none';
axis([-1,1,-1,1,-5,5])
title(['[' num2str(nn) ',' num2str(mm) ']'],'Position',[0.7,1,1],'FontSize',15)
axis off
set(gca,'CameraViewAngle',4)
