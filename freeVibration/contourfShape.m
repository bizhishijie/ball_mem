function s=contourfShape(shape,r,theta,z,z0)
ci=dsearchn(z,z0);
shape_tmp=shape(:,:,ci);
s=contour3(r*cos(theta),r*sin(theta),shape_tmp);
end