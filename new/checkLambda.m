temp=zeros(size(lambda));
for ii=1:size(lambda,1)
    for jj=1:size(lambda,1)
        temp(ii,jj)=sum(lambda(ii,:).*lambda(jj,:));
    end
end