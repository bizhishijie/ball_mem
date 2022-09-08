function OmegaRe=OmegaResonanceNew(n,order_max,cs_ca,L_a,sigma_a_rho, ...
    rootBessel,rootBesselDiff)
% n,order,cs/ca,L/a,sigma/a/rho,gamma*pa/a,root of bessel,root of bessel's derivative
% load('rootBessel.mat')
% load('rootBesselDiff.mat')
% OmegaResonance(0,10, 0.3,1.5, 3.3 ,rootBessel,rootBesselDiff)
% as B in the paper
if n==0% ii从0开始编号
    %lambda
    lambda=zeros(order_max,order_max);
    for ii=1:order_max
        for mm=1:order_max
            lambda(ii,mm)=2*rootBessel(n+1,mm)/(rootBesselDiff(n+1,ii)^2-rootBessel(n+1,mm)^2);
        end
    end

    Omega_loop=0:0.001:50;% as the paper
    det_loop=zeros(1,length(Omega_loop));

    for i_loop=1:length(Omega_loop)
        Omega=Omega_loop(i_loop);
        % A
        A=zeros(order_max);
        for ii=1:order_max%  ii means i,jj means j in the paper
            for jj=1:order_max
                A(ii,jj)=sum(lambda(ii,1:order_max).*lambda(jj,1:order_max)./ ...
                    (cs_ca^2*rootBessel(n+1,1:order_max).^2-Omega^2))/sigma_a_rho;
            end
        end
        % K
        K=zeros(order_max);

        for ii=1:order_max
            sqrt_tmp=sqrt(abs(rootBesselDiff(n+1,ii)^2-Omega^2));
            if Omega<rootBesselDiff(n+1,ii)
                K(ii,ii)=-Omega^2/sqrt_tmp*coth(L_a*sqrt_tmp);% gamma*p_a/a can be eliminate
            elseif Omega>rootBesselDiff(n+1,ii)
                K(ii,ii)=Omega^2/sqrt_tmp*cot(L_a*sqrt_tmp);
            else
                K(ii,ii)=inf;
            end
        end

        mat_temp=K*A;
        %         mat_temp(isnan(mat_temp))=0;
        det_loop(i_loop)=det(eye(order_max)+mat_temp);

    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
else%  ii从1开始编号
    % lambda
    lambda=zeros(order_max,order_max);
    for ii=1:order_max
        for mm=1:order_max
            lambda(ii,mm)=rootBesselDiff(n+1,ii+1)/sqrt(rootBesselDiff(n+1,ii+1)^2-n^2)*...
                2*rootBessel(n+1,mm)/(rootBesselDiff(n+1,ii+1)^2-rootBessel(n+1,mm)^2);
        end
    end

    Omega_loop=0:0.0005:5;% as the paper
    det_loop=zeros(1,length(Omega_loop));

    for i_loop=1:length(Omega_loop)
        Omega=Omega_loop(i_loop);
        % A
        A=zeros(order_max);
        for ii=1:order_max%  ii means i,jj means j in the paper,ii无需+1
            for jj=1:order_max
                A(ii,jj)=sum(lambda(ii,1:order_max).*lambda(jj,1:order_max)./ ...
                    (cs_ca^2*rootBessel(n+1,1:order_max).^2-Omega^2))/sigma_a_rho;
            end
        end
        % K
        K=zeros(order_max);

        for ii=1:order_max
            sqrt_tmp=sqrt(abs(rootBesselDiff(n+1,ii+1)^2-Omega^2));
            if Omega<rootBesselDiff(n+1,ii+1)
                K(ii,ii)=-Omega^2/sqrt_tmp*coth(L_a*sqrt_tmp);% gamma*p_a/a can be eliminate
            elseif Omega>rootBesselDiff(n+1,ii+1)
                K(ii,ii)=Omega^2/sqrt_tmp*cot(L_a*sqrt_tmp);
            else
                K(ii,ii)=inf;
            end
        end
        mat_temp=K*A;
        %         mat_temp(isnan(mat_temp))=0;
        det_loop(i_loop)=det(eye(order_max)+mat_temp);
    end
end

% 取消注释以作图
figure(1);clf;hold on
plot(Omega_loop(det_loop>=0),det_loop(det_loop>=0),'.r')
plot(Omega_loop(det_loop<0),-det_loop(det_loop<0),'.b')
set(gca,'yscale','log')
drawnow

[~,OmegaRe1]=findpeaks(-abs(det_loop),Omega_loop);

tempSign=sign(det_loop(1:end-1)).*sign(det_loop(2:end));
OmegaRe2= Omega_loop(unique([find(tempSign==-1);find(tempSign==-1)+1;find(tempSign==-1)-1]));
OmegaRe=intersect(OmegaRe1,OmegaRe2);
OmegaRe=OmegaRe(1:order_max);

%%%%%%%%%%
trans=zeros(order_max);% 膜腔之间的变换矩阵
% svdMat=zeros(length(OmegaRe),1);
for io =1: order_max
    Omega=OmegaRe(io);
    if n==0
        A=zeros(order_max);
        for ii=1:order_max%  ii means i,jj means j in the paper
            for jj=1:order_max
                A(ii,jj)=sum(lambda(ii,1:order_max).*lambda(jj,1:order_max)./ ...
                    (cs_ca^2*rootBessel(n+1,1:order_max).^2-Omega^2))/sigma_a_rho;
            end%%%%%%%%%%%
        end

        % K
        K=zeros(order_max);

        for ii=1:order_max
            sqrt_tmp=sqrt(abs(rootBesselDiff(n+1,ii)^2-Omega^2));
            if Omega<rootBesselDiff(n+1,ii)
                K(ii,ii)=-Omega^2/sqrt_tmp*coth(L_a*sqrt_tmp);% gamma*p_a/a can be eliminate
            elseif Omega>rootBesselDiff(n+1,ii)
                K(ii,ii)=Omega^2/sqrt_tmp*cot(L_a*sqrt_tmp);
            else
                K(ii,ii)=inf;
            end
        end
    else%%%%%%%%%%%%n>0
        A=zeros(order_max);
        for ii=1:order_max%  ii means i,jj means j in the paper,ii无需+1
            for jj=1:order_max
                A(ii,jj)=sum(lambda(ii,1:order_max).*lambda(jj,1:order_max)./ ...
                    (cs_ca^2*rootBessel(n+1,1:order_max).^2-Omega^2))/sigma_a_rho;
            end
        end
        % K
        K=zeros(order_max);

        for ii=1:order_max
            sqrt_tmp=sqrt(abs(rootBesselDiff(n+1,ii+1)^2-Omega^2));
            if Omega<rootBesselDiff(n+1,ii+1)
                K(ii,ii)=-Omega^2/sqrt_tmp*coth(L_a*sqrt_tmp);% gamma*p_a/a can be eliminate
            elseif Omega>rootBesselDiff(n+1,ii+1)
                K(ii,ii)=Omega^2/sqrt_tmp*cot(L_a*sqrt_tmp);
            else
                K(ii,ii)=inf;
            end
        end
    end
    mat_temp=A*K;
    %     mat_temp(isnan(mat_temp))=0;
    mat_temp=eye(order_max)+mat_temp;
    %     disp(det(mat_temp))
    mat_svd_temp=sort(svd(mat_temp));
    transMat_temp=null(mat_temp,mat_svd_temp(1)+order_max*eps);
    %     svdMat(io)=(mat_svd_temp(2)-mat_svd_temp(1))/mat_svd_temp(1);
    %     if size(transMat_temp,2)~=1
    %         warning('wrong')
    %     end
    trans(io,:)=transMat_temp';
    %      计算并保存变换矩阵
end
trans=trans';% 便于projection中直接乘
save(sprintf('OmegaRe/%d.mat',n),'OmegaRe');
save(sprintf('Trans/%d.mat',n),'trans');
% save(sprintf('svd/%d.mat',n),'svdMat');
end