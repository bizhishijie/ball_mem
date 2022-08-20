function OmegaRe=OmegaResonanceNew(n,order_max,cs_ca,L_a,sigma_a_rho, ...
    rootBessel,rootBesselDiff)
% n，阶数，cs/ca，L/a，sigma/a/rho，gamma*pa/a，bessel函数的根、极值点
% load('rootBessel.mat')
% load('rootBesselDiff.mat')
% OmegaResonanceNew(0,50, 0.6,5, 3.3 ,rootBessel,rootBesselDiff)
lambda=zeros(order_max,order_max);
% n的最低为0，只有在rootBessel与rootBesselDiff中需要+1
if n==0% ii对应i，mm对应m
    for ii=1:order_max
        for mm=1:order_max
            lambda(ii,mm)=2*rootBessel(n+1,mm)/(rootBesselDiff(n+1,ii)^2-rootBessel(n+1,mm)^2);
        end
    end
else
    for ii=1:order_max
        for mm=1:order_max
            lambda(ii,mm)=rootBesselDiff(n+1,mm)/sqrt(rootBesselDiff(n+1,mm)^2-n^2)*...
                2*rootBessel(n+1,mm)/(rootBesselDiff(n+1,ii)^2-rootBessel(n+1,mm)^2);
        end
    end
end
lambda(isnan(lambda)) = 0;% lambda近似满足正交归一性，见checkLambda
% Omega_loop=linspace(0,30,10000)
Omega_loop=0:0.001:10;% 见论文
det_loop=zeros(1,length(Omega_loop));

for i_loop=1:length(Omega_loop)

    Omega=Omega_loop(i_loop);
    A=zeros(order_max,order_max);
    for ii=1:order_max% ii对应i，jj对应j，与上不同，对应论文中
        for jj=1:order_max
            A(ii,jj)=sum(lambda(ii,:).*lambda(jj,:)./ ...
                (cs_ca^2*rootBessel(n+1,1:order_max).^2-Omega^2));
        end
    end
    A=A/sigma_a_rho;

    K=zeros(order_max,order_max);
    for ii=1:order_max
        sqrt_tmp=sqrt(abs(rootBesselDiff(n+1,ii)^2-Omega^2));
        if Omega<rootBesselDiff(n+1,ii)
            K(ii,ii)=-Omega^2/sqrt_tmp*coth(L_a*sqrt_tmp);% gamma*p_a/a 略去
        elseif Omega>rootBesselDiff(n+1,ii)
            K(ii,ii)=Omega^2/sqrt_tmp*cot(L_a*sqrt_tmp);
        else
            K(ii,ii)=inf;
        end
    end

    det_loop(i_loop)=det(1+K*A);
end

figure(1);clf;hold on
plot(Omega_loop(det_loop>=0),det_loop(det_loop>=0),'.r')
plot(Omega_loop(det_loop<0),-det_loop(det_loop<0),'.b')
set(gca,'yscale','log')

tempSign=sign(det_loop(1:end-1)).*sign(det_loop(2:end));
OmegaRe=Omega_loop(tempSign==-1);
save(sprintf('Omega/%d.mat',n),'OmegaRe');
end