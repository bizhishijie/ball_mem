function gui_draw_all
shape_c=cell(order_max+1,order_max);
for nn=0:order_max
    for mm=1:order_max
        if nn==0 && mm==1
            z1=sqrt(2)*ones(size(r))*cos(nn*theta)/sqrt(pi)/sqrt(2);
        elseif nn==0
            z1=sqrt(2)*rootBesselDiff(nn+1,mm)*besselj(nn,rootBesselDiff(nn+1,mm)*r)/...
                besselj(nn,rootBesselDiff(nn+1,mm))/sqrt(rootBesselDiff(nn+1,mm)^2-nn^2)*...
                cos(nn*theta)/sqrt(pi)/sqrt(2);
        else%mm 需要+1
            z1=sqrt(2)*rootBesselDiff(nn+1,mm+1)*besselj(nn,rootBesselDiff(nn+1,mm+1)*r)/...
                besselj(nn,rootBesselDiff(nn+1,mm+1))/sqrt(rootBesselDiff(nn+1,mm+1)^2-nn^2)*...
                cos(nn*theta)/sqrt(pi);
        end
        % 先是r,后是theta
        shape_c{nn+1,mm}=z1;
    end
end

shape_m=cell(order_max+1,order_max);
for nn=0:order_max
    for mm=1:order_max
        z1=sqrt(2)*besselj(nn,rootBessel(nn+1,mm)*r)/...
            besselj(nn+1,rootBessel(nn+1,mm))*...
            cos(nn*theta)/sqrt(pi)/sqrt(1+(nn==0));
        % 先是r,后是theta
        shape_m{nn+1,mm}=z1;
    end
end

a=1;ca=340;
omega=Omega*ca/a;
deep=5;

cnt=0;
mkdir(sprintf('./couple/fig/%d.%d/',n,m));
z_max=max(max(abs(shape_c{n+1,m})));
omega_tmp=omega(n+1,m);

for t=linspace(0,2*pi/omega_tmp,100)
    cnt=cnt+1;
    s1=surfShape(app.UIAxes2,shape_c{n+1,m}*cos(omega_tmp*t)-deep,r,theta);
    s1.EdgeColor="none";
    app.UIAxes2.NextPlot='add';
    s2=surfShape(app.UIAxes2,shape_m{n+1,m}*cos(omega_tmp*t),r,theta);
    s2.EdgeColor="none";
    axis(app.UIAxes2,[-1,1,-1,1,-20,20]);
    load('s_w0.mat');
    load('s_h0.mat');
    view(app.UIAxes2,see_where0,see_high0);  % 设置视点位置
    axis off
    drawnow
    pause(0.05)
    app.UIAxes2.NextPlot='replacechildren';
    %saveas(gcf,sprintf("./couple/fig/%d.%d/%d.bmp",n,m,cnt))

end
end