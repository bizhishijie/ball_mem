clf;close;
figure(1)
% uiopen('G:\p14\7.17实验\不同半径&落点\频谱\高度4cm，白色薄膜0.095\A2_0.5 频谱.fig',1)
% h1=findall(gca,'type','line');
% x1=get(h1,'xdata');
% y1=get(h1,'ydata');
uiopen('G:\p14\7.17实验\不同半径&落点\频谱\高度4cm，白色薄膜0.095\A2_0.75 频谱.fig',1)
h2=findall(gca,'type','line');
x2=get(h2,'xdata');
y2=get(h2,'ydata');
% figure(2)
% plot(x1,y1,'LineWidth',1)
% axis([0 1200 0 1])
% hold on
plot(x2,y2,'LineWidth',1)
axis([0 1200 0 1])
hold on
m=4;
for iiiii=1:m
load(['F:\final\OmegaRe\' num2str(iiiii-1) '.mat'])
F=linspace(1,3000,3000);
A=OmegaRe*ca/a/2/pi;
A1=floor(A);
A2=zeros(1,3000);
A2(A1)=1;
figure(2)
plot(F,A2(1:3000),'Color',[0 1-iiiii/m 1])
hold on
axis([0 1200 0 1])
end
