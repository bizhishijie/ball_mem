function [f,w]=efft1(varargin)
% efft1(X,Y,type,number)被用于求解离散数值的一维傅里叶变换
% X,Y是具有同样长度的一维矩阵,其中X不要求等间距
% [xx,yy]=efft1(X,Y) 求解由X,Y组成的傅里叶变换，即把Y沿着X进行傅里叶变换
% [xx,yy]=efft1(X,Y,type,number)
% 求解由X,Y组成的傅里叶变换，即把Y沿着X进行傅里叶变换,同时指定傅里叶变换结果是否需要归一化、傅里叶变换的基础数据
% type type的可选择值有'0','-1',缺省三种，'0'表示结果不需要归一化，'-1'表示结果需要归一化,缺省表达时，默认归一化
% number number的可选择值是任意以字符串表示的数字，例如'10000'，或缺省两种
% 缺省表达时，默认为X的全部长度
% 当提示'more number points were needed'时，请适当增加输入的number数值，或者选择缺省表达
% 当提示'the sampling points number were banned over the length of the data
% that you import'时，请适当减小输入的number数值，或者选择缺省表达
% 最终输出的横轴已经自动修改为频率

% 杨敬整理、夏成杰修改、宫殿锦丰提供原始函数
args = matlab.images.internal.stringToChar(varargin);
[X,Y,type,number] = parse_inputs(args{:});
if strcmp(number,'all')
    numb=length(X);
else
    numb=str2num(number);
end
if length(find(abs(diff(diff(X)))<10^(-5)))==length(X)-2
    
else
    [X,XI]=unique(X);
    Yi=Y(XI);
    Xi=linspace(min(X),max(X),numb);
    Yi=interp1(X,Yi,Xi,'linear');
    X=Xi;
    Y=Yi;
end
if length(find(abs(diff(diff(X)))<10^(-5)))==length(X)-2                          %此时横坐标等间距
    deltax=diff(X);
    f=linspace(0,1/(deltax(1)),numb);
    if str2num(type)==-1
        w=abs(fft(Y,numb)/sum(abs(Y(1:numb))));
    else
        w=abs(fft(Y,numb));
    end
end
end

function [X,Y,type,number] =parse_inputs(varargin)
X=varargin{1};
Y=varargin{2};
Type={'0','-1'};
type='-1';
if nargin>=3
    for ii=3:nargin
        if ischar(varargin{ii})
            str = lower(varargin{ii});
            j = find(strcmp(str,Type));
            k = str2num(str);
            if ~isempty(j) && (k==0 || k==-1)
                type=Type{j(1)};
                number='all';
                continue
            else
                if k>=20 && k<length(Y)
                    number=num2str(k);
                elseif k>=length(Y)
                    error('the sampling points number were banned over the length of the data that you import');
                    
                else
                    error('more number points were needed');
                end
            end
        end
    end
else
    type='-1';
    number='all';
end
end