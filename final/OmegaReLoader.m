fileList=dir('./OmegaRe/*.mat');
order_max=10;% 需要修改%%%%%%%%%%%%%%%
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

% order_max=5;% 需要修改
transMat=zeros(length(fileList),order_max,order_max);
for ni =1:length(fileList)
    load(strcat(fileList(ni).folder,'/', fileList(ni).name));
    transMat(ni,:,:)=trans(1:order_max,1:order_max);
    % 报错就增大Omega_loop
%     det(trans(ni));% 这个数字很接近0，计算应该是对的
end
save('transMat.mat',"transMat")