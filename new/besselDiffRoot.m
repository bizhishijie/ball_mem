function rootBesselDiff  =besselDiffRoot(n,N) 
% besselDiffRoot(50,50) 
% 贝塞尔函数的导数的根
r=0:0.001:1e3;
% besTemp=zeros(size(r));
rootBesselDiff = zeros(n+1,N);
for ni=0:n
    besTemp=besselj(ni,r);
    [~,locs]=findpeaks(abs(besTemp));
    if length(locs)>N
        locs=locs(1:N);
    end
    rootBesselDiff(ni+1,:)=r(locs);
end
rootBesselDiff(1,:)=[0,rootBesselDiff(1,1:end-1)];
save("rootBesselDiff.mat",'rootBesselDiff')
end