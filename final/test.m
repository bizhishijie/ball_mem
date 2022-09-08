% a=0.05;
% L=0.1;
% rho_mem=0.095;
% rho_air=1.184;
% a0=0.1*a;
% T=41.9;
% cs=sqrt(T/rho_mem);
% ca=345;
% 
% sigma = 0.120;

% cs_ca=cs/ca;L_a=L/a;sigma_a_rho=sigma/a/rho_air;

load('rootBessel.mat')
load('rootBesselDiff.mat')
load('OmegaRe.mat')

a=0.047; %半径
L = 0.1;% 鼓的深度
gamma = 1.4; % 空气绝热系数
T=19; % 张力19.9
rho = 1.184;% 空气密度
ca = 345;% 空气中声速
sigma = 0.120;% 气球面密度
pa = 1e5;% 大气压强
cs = sqrt(T/sigma);% 橡胶中波的传播速度
cs_ca=cs/ca;L_a=L/a;sigma_a_rho=sigma/a/rho;

for n=0:20
    disp(n)
    OmegaResonanceNew(n,20, cs_ca,L_a, sigma_a_rho,rootBessel,rootBesselDiff);
end

fileList=dir('./OmegaRe/*.mat');
order_max=20;% 需要修改%%%%%%%%%%%%%%%
Omega=zeros(length(fileList),order_max);
for ni =1:length(fileList)
    load(strcat(fileList(ni).folder,'/', fileList(ni).name));
    while length(OmegaRe)<order_max
        OmegaRe=[OmegaRe nan];
    end 
    Omega(ni,:)=OmegaRe(1:order_max);
end
save('OmegaRe.mat',"Omega")

fileList=dir('./Trans/*.mat');

transMat=zeros(length(fileList),order_max,order_max);
for ni =1:length(fileList)
    load(strcat(fileList(ni).folder,'/', fileList(ni).name));
    transMat(ni,:,:)=trans(1:order_max,1:order_max);
    % 报错就增大Omega_loop
%     det(trans(ni));% 这个数字很接近0，计算应该是对的
end
save('transMat.mat',"transMat")

n=1;%%%%%%%%%%%%%% 请改这个
% rootBesselDiff=rootBesselDiff(n,:);
% rootBessel=rootBessel(1,:);
Omega=Omega(n,:);
f_cav_close=sqrt(rootBesselDiff(n,:).^2+pi^2/L_a^2);
f_cav_open=sqrt(rootBesselDiff(n,:).^2+(pi/2)^2/L_a^2);
f_mem=rootBessel(n,:)*cs_ca;
f_cou=Omega;

f_cav_close=sort(f_cav_close(:));
f_cav_open=sort(f_cav_open(:));
f_mem=sort(f_mem(:));
f_cou=sort(f_cou(:));

clf;hold on

plot(f_cav_close,4*ones(size(f_cav_close)),'k*')
plot(f_cav_open,3*ones(size(f_cav_open)),'g*')
plot(f_mem,2*ones(size(f_mem)),'r*')
plot(f_cou,ones(size(f_cou)),'b*')

axis([0.3 10 0 5]);
grid on
set(gca,'xscale','log')