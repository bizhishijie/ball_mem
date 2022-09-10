%https://ww2.mathworks.cn/help/matlab/ref/plot.html?s_tid=srchtitle
%�õط���
%%
x=linspace(0,2*pi,1000);
y=sin(x);
plot(x,y)%��ǰһ������Ϊ�����꣬��һ������Ϊ�����꣬���㡣
%x��y����ӵ����ͬ��Ԫ������ͬΪ��������������
%%
clf%�����ǰͼ��
%%
x=linspace(-2*pi,2*pi);
y1 = sin(x);
y2 = cos(x);
plot(x,y1,x,y2)%������x��y��ͬʱ���ƶ�������
%%
close%�رյ�ǰ����
%%
figure%����һ���´���
%%
axis([0 2*pi -1 1])%���û�ͼ�������᷶Χ����[]����
%ǰ������Ϊ�����귶Χ��������Ϊ�����귶Χ
%range=[0 2*pi -1 1]%Ҳ�ɽ�[]�ڵķ�Χ������Ϊ������Ч����ͬ
%axis(range)
%%
axis equal%����������ı�����Ϊ1��1

%%
box on%��ʾ�ң���������
%%
grid on%��ʾ������
%%
clf
x = 0:pi/100:2*pi;
y1 = sin(x);
y2 = sin(x-0.25);
y3 = sin(x-0.5);
y4 = sin(x-0.75);
plot(x,y1,'ro')%��ÿ��xy����á���������һ�����ߵĸ�ʽ,
%r�����������Ϊ��ɫ��o������Ȧ��ʾÿ����
hold on%ִ���µ�plot����ʱ�����ԭ�ȵ�����
plot(x,y2,'b.')%bΪ��ɫ���Դ����ƣ�.Ϊ���Ե���ʾÿ����
plot(x,y3,'g*')%*Ϊ��*��ʾÿ����
plot(x,y4,'--o')%--�����������,��һ�λ��ƶ��ַ���
%%
clf
plot(x,y1,'o')%
hold on%
plot(x,y2,'.')%
plot(x,y3,'*')%
plot(x,y4,'--o')%
%���https://ww2.mathworks.cn/help/matlab/ref/plot.html?s_tid=srchtitle#btzitot-LineSpec
%������ÿ�����ߵ���ɫ����matlab���õ���ɫ��
%��plot������Ⱥ�˳���趨ÿ�����ߵ���ɫ
lg=legend('y1','y2','y3','y4')%�Դ�����ÿ���ߵ�ͼ��,lg=Ϊ���þ�����������
set(lg,'box','off')%ȥ��ͼ�������ᣬ�Ƽ�
%���https://ww2.mathworks.cn/help/matlab/ref/legend.html?searchHighlight=legend&s_tid=srchtitle
%%
x = linspace(0,10);
y = sin(x);
stp=1:5:length(y);
plot(x,y,'-o','MarkerIndices',stp)
%��ÿ�����ߺ��ԡ�����עҪ���ĵı�������������ű���ֵ
%'MarkerIndices'Ϊֻ��ʾ��stp�ĵ�
%%
clf
x = -pi:pi/10:pi;
y = tan(sin(x)) - sin(tan(x));
plot(x,y,'--gs',...%g��ʾ��ɫ��s��ʾ�Է������ÿ���㣬������o��.��ͬ
    'LineWidth',2,...%�����߿�
    'MarkerSize',10,...%���õ�Ĵ�С
    'MarkerFaceColor',[0.5,0.5,0.5],...%���õ��ڲ�����ɫ���Թ�һ����[R,G,B]����
    'MarkerEdgeColor','b')%���õ��Ե����ɫ,bΪ��ɫ����ɫΪk
%����ֵ����ĸ����ɫ���ö����Ա�[R,G,B]����
%��xy������ġ�--gs���ڲ�����ɫ������[R,G,B]����
%...��ʾ��������δ������ת����һ��
%%
clf
x = linspace(0,10,500);
y = cos(5*x);
plot(x,y,'Color',[0,0.7,0.9],'LineWidth',1.5)%�á�color�������������ߵ���ɫ
title('2-D Line Plot')%���ñ���
xlabel('x')%����x������
ylabel('cos(5x)')%����y������
%%
set(gca,'LineWidth',2,'FontName','΢���ź�','FontSize',20)
%      �������߿�     ����������        �������ֺ�
title('2-D Line Plot','FontName','΢���ź�','FontSize',28)
%fontname�ȶ����ֵ����ö���title��xlable��ͬ����Ч
%�����������õĻ�title������������̶�������ͬ
%title��ʹ�����https://ww2.mathworks.cn/help/matlab/ref/title.html?searchHighlight=title&s_tid=srchtitle
%%
clf
x=1:0.1:1e5;
y=x.^2;
plot(x,y)
%%
set(gca,'XScale','log','YScale','log')%��x��y���Ϊ��������
%��log����Ϊlinear���ԸĻ���������
grid on
%% ���±༭�ѻ��Ƶ����ߣ�ʹ�þ�����Ƽ�ʹ��
clf
x = linspace(-2*pi,2*pi);
y1 = sin(x);
y2 = cos(x);
p = plot(x,y1,x,y2);%p�д������������ߵ���Ϣ�����зֱ���б༭
%%
p(1).LineWidth = 2;%���÷�������þ���Ԫ����ͬ
%��ʵ��p����һ��Ԫ��Ϊֱ�ߵľ�����������������ʹ�ýṹ��ķ�����ͬ
p(2).Marker = '*';
%%
delete(p(1))%������delete����ɾ���Ѿ����õ�����
%%
set(gcf,'Position',[0 0 800 500])%����ͼƬ��ʾ��λ�úʹ�С
%[]�ڷֱ�Ϊx����㡢y����㡢x���յ㡢y���յ�
%% ��ά���ߵĻ��ƣ�ǰ�����й��ܾ���ʹ��
x=linspace(0,2*pi,1000);
y=linspace(0,2*pi,1000);
z=sin(x+y);
plot3(x,y,z,'-o','MarkerIndices',1:100:length(z),'LineWidth',2,'MarkerSize',10)
grid on
box on
set(gca,'LineWidth',2,'FontName','΢���ź�','FontSize',20)
set(gcf,'Position',[0 0 1000 800])
axis equal
xlabel('x')
ylabel('y')
zlabel('z')
title('��άͼ��')
view(45,30)%���ù۲�ĽǶȣ�(ʱ��,����)���ĸ����ԾͶ��ˣ�������0��45
%% �������ܸ�ɶ��
clf
f=@(t,Y)[10*(-Y(1)+Y(2));...
28*Y(1)-Y(2)-Y(1)*Y(3);...
Y(1)*Y(2)-8*Y(3)/3;];
[t,Y]=ode45(f,[0,200],[12,2,9]);
plot3(Y(:,1),Y(:,2),Y(:,3),'color',[255/255,200/255,60/255])
set(gca,'color','k')%��������������ɫ
set(gcf,'color','k')%��������������ɫ
grid on
box on
set(gca,'GridColor','w')%����������ɫ
set(gca,'xcolor','w','ycolor','w','zcolor','w');%������������ɫ
set(gca,'LineWidth',2,'FontName','΢���ź�','FontSize',20)
set(gcf,'Position',[0 0 1000 800])
title('������������','color','w','fontsize',30)
view(30,10)