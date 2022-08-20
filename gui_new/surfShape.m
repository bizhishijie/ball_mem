function s=surfShape(fig,shape,r,theta)
a=length(r);b=length(theta);
theta=[theta theta(1)];
shape_tmp=zeros(a,b+1);
shape_tmp(1:a,1:b)=shape;
% shape_tmp(end,:)=shape_tmp(1,:);
shape_tmp(:,end)=shape_tmp(:,1);
s=surf(fig,r*cos(theta),r*sin(theta),shape_tmp);
end