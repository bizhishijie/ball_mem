%https://ww2.mathworks.cn/help/matlab/ref/plot.html?s_tid=srchtitle
%好地方↑
%%
x=linspace(0,2*pi,1000);
y=sin(x);
plot(x,y)%以前一个变量为横坐标，后一个变量为纵坐标，画点。
%x和y必须拥有相同的元素数且同为行向量或列向量
%%
clf%清除当前图层
%%
x=linspace(-2*pi,2*pi);
y1 = sin(x);
y2 = cos(x);
plot(x,y1,x,y2)%输入多个x，y可同时绘制多条曲线
%%
close%关闭当前窗口
%%
figure%创建一个新窗口
%%
axis([0 2*pi -1 1])%设置画图的坐标轴范围，用[]框起
%前两个数为横坐标范围，后两个为纵坐标范围
%range=[0 2*pi -1 1]%也可将[]内的范围单独设为变量，效果相同
%axis(range)
%%
axis equal%将横纵坐标的比例设为1：1

%%
box on%显示右，上坐标轴
%%
grid on%显示网格线
%%
clf
x = 0:pi/100:2*pi;
y1 = sin(x);
y2 = sin(x-0.25);
y3 = sin(x-0.5);
y4 = sin(x-0.75);
plot(x,y1,'ro')%在每组xy后可用‘’设置这一条曲线的格式,
%r代表词条曲线为红色，o代表用圈显示每个点
hold on%执行新的plot命令时不清楚原先的曲线
plot(x,y2,'b.')%b为蓝色，以此类推；.为将以点显示每个点
plot(x,y3,'g*')%*为以*显示每个点
plot(x,y4,'--o')%--代表绘制虚线,可一次绘制多种符号
%%
clf
plot(x,y1,'o')%
hold on%
plot(x,y2,'.')%
plot(x,y3,'*')%
plot(x,y4,'--o')%
%详见https://ww2.mathworks.cn/help/matlab/ref/plot.html?s_tid=srchtitle#btzitot-LineSpec
%不设置每条曲线的颜色是以matlab内置的颜色库
%按plot命令的先后顺序设定每条曲线的颜色
lg=legend('y1','y2','y3','y4')%以此设置每条线的图例,lg=为设置句柄，详见后文
set(lg,'box','off')%去掉图例坐标轴，推荐
%详见https://ww2.mathworks.cn/help/matlab/ref/legend.html?searchHighlight=legend&s_tid=srchtitle
%%
x = linspace(0,10);
y = sin(x);
stp=1:5:length(y);
plot(x,y,'-o','MarkerIndices',stp)
%在每条曲线后以‘’标注要更改的变量，后面紧跟着变量值
%'MarkerIndices'为只显示第stp的点
%%
clf
x = -pi:pi/10:pi;
y = tan(sin(x)) - sin(tan(x));
plot(x,y,'--gs',...%g表示绿色，s表示以方块绘制每个点，作用与o，.相同
    'LineWidth',2,...%设置线宽
    'MarkerSize',10,...%设置点的大小
    'MarkerFaceColor',[0.5,0.5,0.5],...%设置点内部的颜色，以归一化的[R,G,B]设置
    'MarkerEdgeColor','b')%设置点边缘的颜色,b为蓝色，黑色为k
%变量值的字母类颜色设置都可以被[R,G,B]代替
%但xy后紧跟的‘--gs’内部的颜色不能用[R,G,B]代替
%...表示此行命令未结束，转接下一行
%%
clf
x = linspace(0,10,500);
y = cos(5*x);
plot(x,y,'Color',[0,0.7,0.9],'LineWidth',1.5)%用‘color’来自由设置线的颜色
title('2-D Line Plot')%设置标题
xlabel('x')%设置x轴名称
ylabel('cos(5x)')%设置y轴名称
%%
set(gca,'LineWidth',2,'FontName','微软雅黑','FontSize',20)
%      坐标轴线宽     坐标轴字体        坐标轴字号
title('2-D Line Plot','FontName','微软雅黑','FontSize',28)
%fontname等对文字的设置对于title，xlable等同样有效
%若不单独设置的话title等体与坐标轴刻度字体相同
%title的使用详见https://ww2.mathworks.cn/help/matlab/ref/title.html?searchHighlight=title&s_tid=srchtitle
%%
clf
x=1:0.1:1e5;
y=x.^2;
plot(x,y)
%%
set(gca,'XScale','log','YScale','log')%将x、y轴改为对数坐标
%将log设置为linear可以改回先行坐标
grid on
%% 重新编辑已绘制的曲线（使用句柄）推荐使用
clf
x = linspace(-2*pi,2*pi);
y1 = sin(x);
y2 = cos(x);
p = plot(x,y1,x,y2);%p中储存了两条曲线的信息，可有分别进行编辑
%%
p(1).LineWidth = 2;%调用方法与调用矩阵元素相同
%事实上p就是一个元素为直线的矩阵，这与其他语言中使用结构体的方法相同
p(2).Marker = '*';
%%
delete(p(1))%可以用delete命令删除已经画好的曲线
%%
set(gcf,'Position',[0 0 800 500])%调整图片显示的位置和大小
%[]内分别为x轴起点、y轴起点、x轴终点、y轴终点
%% 三维曲线的绘制，前述所有功能均可使用
x=linspace(0,2*pi,1000);
y=linspace(0,2*pi,1000);
z=sin(x+y);
plot3(x,y,z,'-o','MarkerIndices',1:100:length(z),'LineWidth',2,'MarkerSize',10)
grid on
box on
set(gca,'LineWidth',2,'FontName','微软雅黑','FontSize',20)
set(gcf,'Position',[0 0 1000 800])
axis equal
xlabel('x')
ylabel('y')
zlabel('z')
title('三维图像')
view(45,30)%设置观察的角度，(时角,俯角)，改改试试就懂了，来几个0和45
%% 这玩意能干啥？
clf
f=@(t,Y)[10*(-Y(1)+Y(2));...
28*Y(1)-Y(2)-Y(1)*Y(3);...
Y(1)*Y(2)-8*Y(3)/3;];
[t,Y]=ode45(f,[0,200],[12,2,9]);
plot3(Y(:,1),Y(:,2),Y(:,3),'color',[255/255,200/255,60/255])
set(gca,'color','k')%更改坐标轴内颜色
set(gcf,'color','k')%更改坐标轴外颜色
grid on
box on
set(gca,'GridColor','w')%更改网格颜色
set(gca,'xcolor','w','ycolor','w','zcolor','w');%更改坐标轴颜色
set(gca,'LineWidth',2,'FontName','微软雅黑','FontSize',20)
set(gcf,'Position',[0 0 1000 800])
title('洛伦兹吸引子','color','w','fontsize',30)
view(30,10)