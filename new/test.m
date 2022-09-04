open('A1-4.7cm-距膜圆心0mm 频谱.fig')
load('./OmegaRe/0.mat')
h=zeros(1,13000);
ca=346;% 空气中声速
a = 0.047;% 鼓的半径
h(floor(OmegaRe*ca/a/2/pi))=1;  
hold on
plot(h)
axis([0 2000 0 1])