function rootBessel  =besselRoot(n,N) 
% besselRoot(50,51) 
% 贝塞尔函数的导数的根
r=0:0.00001:250;
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
save("rootBessel.mat",'rootBessel')
end