function rootBessel  =besselRootNew(n,N) 
% besselRootNew(50,50) 
% 贝塞尔函数的导数的根
r=0:0.001:1e3;
% besTemp=zeros(size(r));
rootBessel = zeros(n+1,N);
for ni=0:n
    besTemp=besselj(ni,r);
    [~,locs]=findpeaks(-abs(besTemp));
    if length(locs)>N
        locs=locs(1:N);
    end
    rootBessel(ni+1,:)=r(locs);
end
rootBessel(1,:)=[0,rootBessel(1,1:end-1)];
save("rootBessel.mat",'rootBessel')
end