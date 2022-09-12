for nn=1:3
    for mm=1:3
        filePath=['./pic/' num2str(nn) '.' num2str(mm)];
        fileList=dir([filePath '/*.jpg']);
        for i = 1:length(fileList)-1
            str = [filePath '/' num2str(i) '.jpg'];
            disp(str)
            A = imread(str);
            [I, map] = rgb2ind(A, 256);

            if (i == 1)
                imwrite(I, map,[num2str(nn) '.' num2str(mm) '.gif'], 'DelayTime', 0.1, 'LoopCount', Inf)
            else
                imwrite(I, map,[num2str(nn) '.' num2str(mm) '.gif'], 'WriteMode', 'append', 'DelayTime', 0.1)
            end
        end
    end
end