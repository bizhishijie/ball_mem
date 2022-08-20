function OmegaRe=OmegaResonance(n,order_max,cs_ca,L_a,sigma_a_rho,gamma_p_a,rootBessel,rootBesselDiff)
% n，阶数，cs/ca，深度比半径，sigma/a/rho，gamma*pa/a，bessel函数的根、极值点
% load('rootBessel.mat')
% load('rootBesselDiff.mat')
% OmegaResonance(1,40,5.29,1, 0.0025,700000,rootBessel,rootBesselDiff)
lambda=zeros(order_max,order_max);
for ii=1:order_max
    for jj=1:order_max
        lambda(ii,jj)=2*rootBessel(n+1,jj)/(rootBesselDiff(n+1,ii)^2-rootBessel(n+1,jj)^2);
    end
end

% Omega_loop=linspace(0,30,10000)
Omega_loop=0:0.001:10;
det_loop=zeros(1,length(Omega_loop));

for i_loop=1:length(Omega_loop)
    
    Omega=Omega_loop(i_loop);
    A=zeros(order_max,order_max);
    for ii=1:order_max
        for jj=1:order_max
            A(ii,jj)=sum(lambda(ii,:).*lambda(jj,:)./(cs_ca^2*rootBessel(n+1,1:order_max).^2-Omega^2));
        end
    end
    A=A/gamma_p_a/sigma_a_rho;
    
    K=zeros(order_max,order_max);
    for ii=1:order_max
        sqrt_tmp=sqrt(abs(rootBesselDiff(n+1,ii)^2-Omega^2));
        if Omega<rootBesselDiff(n+1,ii)
            K(ii,ii)=-gamma_p_a*Omega^2/sqrt_tmp*coth(L_a*sqrt_tmp);
        elseif Omega==rootBesselDiff(n+1,ii)
            K(ii,ii)=inf;
        else
            K(ii,ii)=gamma_p_a*Omega^2/sqrt_tmp*cot(L_a*sqrt_tmp);
        end
    end
    
    det_loop(i_loop)=det(1+K*A);
end

figure(1);clf;hold on
plot(Omega_loop(det_loop>=0),det_loop(det_loop>=0),'.r')
plot(Omega_loop(det_loop<=0),-det_loop(det_loop<=0),'.b')
set(gca,'yscale','log')

tempSign=sign(det_loop(1:end-1)).*sign(det_loop(2:end));
OmegaRe=Omega_loop(tempSign==-1);
save(sprintf('Omega/%d.mat',n),'OmegaRe');
end