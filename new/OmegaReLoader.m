fileList=dir('./Omega/*.mat');
order_max=50;
Omega=zeros(length(fileList),order_max);
for ni =1:length(fileList)
    load(strcat(fileList(ni).folder,'/', fileList(ni).name));
    Omega(ni,:)=OmegaRe(1:order_max);
end
save('OmegaRe.mat',"Omega")