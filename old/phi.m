function p=phi(n,m,r,rootBessel)
mu_n_m=rootBessel(n,m);% 贝塞尔方程的根的前数部份
mu1=repmat(mu_n_m,1,1,length(r));%将其第三个维度重整为r的形式
mur1=mu1.*r;
p=zeros(size(mur1));
for ni=n
    p(ni,:,:)=besselj(ni-1,mur1(ni,:,:))./besselj(ni-1,mu1(ni,:,:));
end
end