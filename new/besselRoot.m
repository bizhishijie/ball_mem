function rootBessel  =besselRoot(n,N)    
% 已被替代
% besselRoot(50,50)    
% 求解n阶贝塞尔函数的零点(0-9)
% n为贝塞尔函数阶数 
% N为要求的零点数目
% 每一行为每一阶的若干零点
j = zeros(n+1, N);    % 贝塞尔函数的根
incr = 4.0;
for v = 0 : n
   h = v + 1.9*v^(1/3)+1;
   if (v == 0)             % 0阶贝塞尔函数的第一个零点
       j(v+1,1) = fzero(@(x)besselj(v,x),2);
   else                    % 1阶及以上阶贝塞尔函数的第一个零点
       j(v+1,1) = fzero(@(x)besselj(v,x),h);
   end
   for s = 2 : N           % 贝塞尔函数的第2个及后面的零点
       j(v+1,s) = fzero(@(x)besselj(v,x),j(v+1,s-1)+incr);
   end    
end

rootBessel = j; 

save("rootBessel.mat",'rootBessel')
end