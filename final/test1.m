clf
[y,Fs] = audioread('.\aac\新录音 76.m4a');
y=y(3.3e4:4.3e4,1);
x=1:length(y);
x=x/Fs;
[f,w]=efft1(x,y);
plot(f,w)
xlim([0 1000])
%%
figure
[y,Fs] = audioread('.\aac\新录音 76.m4a');
y=y(4.3e4:6e4,1);
x=1:length(y);
x=x/Fs;
[f,w]=efft1(x,y);
plot(f,w)
xlim([0 1000])
%%
load('.\OmegaRe\0.mat')
OmegaRe0=OmegaRe;
load('.\OmegaRe\1.mat')
OmegaRe1=OmegaRe;
load('.\OmegaRe\2.mat')
OmegaRe2=OmegaRe;
load('.\OmegaRe\3.mat')
OmegaRe3=OmegaRe;
OmegaRe=[ OmegaRe0 OmegaRe1];
OmegaRe=OmegaRe*cs/a;
% 对的波速大约是27
A=zeros(size(y));
A(ceil(OmegaRe))=1;
hold on 
plot(A)

xlim([0 1000])