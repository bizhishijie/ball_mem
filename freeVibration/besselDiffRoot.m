function rootBesselDiff  =besselDiffRoot(n,N)
% besselDiffRoot(50,51)
% 贝塞尔函数的导数的根
r=0:0.00001:250;
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
for ii=2:n+1
    rootBesselDiff(ii,:)=[nan,rootBesselDiff(ii,1:end-1)];
end
save("rootBesselDiff.mat",'rootBesselDiff')
end