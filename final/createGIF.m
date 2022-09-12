fileList=dir('./pic/*.jpg');
for i = 1:length(fileList)-1
    str = ['./pic/' num2str(i) '.jpg'];
    disp(str)
    A = imread(str);
    [I, map] = rgb2ind(A, 256);

    if (i == 1)
        imwrite(I, map, 'movefig.gif', 'DelayTime', 0.1, 'LoopCount', Inf)
    else
        imwrite(I, map, 'movefig.gif', 'WriteMode', 'append', 'DelayTime', 0.1)
    end

end
