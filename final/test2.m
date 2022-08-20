clf
[y,Fs] = audioread('.\aac\8月5日 18点16分.aac');
y=y(2.05e5:2.25e5,1);
x=1:length(y);
x=x/Fs;
[f,w]=efft1(x,y);
plot(f,w)
xlim([0 1000])
hold on
%%
load('.\OmegaRe\0.mat')
OmegaRe0=OmegaRe;
load('.\OmegaRe\1.mat')
OmegaRe1=OmegaRe;
load('.\OmegaRe\2.mat')
OmegaRe2=OmegaRe;
OmegaRe=[ OmegaRe0 OmegaRe1 OmegaRe2];
OmegaRe=OmegaRe*520;
% 对的波速大约是27
a=zeros(size(y));
a(ceil(OmegaRe))=1;
% hold on 
plot(a)

xlim([0 1000])