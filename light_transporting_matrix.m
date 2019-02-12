c=299792458; %光速
epsilon0=8.854187817*power(10,-12); %真空介电常数ε0 = 8.854187817*10-12F/m
E0in=1;   %定义入射到玻璃、第一薄层界面的光电场的大小，默认设为一
layers=4;  %定义器件共有多少个薄层
epsilon=cell(layer,1);  %创建用于存放器件各层介电函数（复折射率）的元胞矩阵
glass=load('glass.txt');  %读取玻璃的介电函数
for j=1:layer
    tmp=[num2str(j),'_epsilon.txt']; %读取各层的介电函数（复折射率）
    epsilon{j}=load(tmp);
end
thickness=load('thickness.txt');  %读入各层的厚度值
inmode=1; %入射模式，1表示垂直入射，2表示偏角p偏振，3表示偏角s偏振
I=cell(size(glass,1),layer+1);
L=cell(size(glass,1),layer); 
Zetaelec=zeros(size(glass,1),layer);
%分别创建空的界面矩阵组和层矩阵组，以及各层的zeta函数值待放数据
%然后对整个器件的光入射方向的坐标轴划分分度值，即取好x序列
interval=0.5;  %设定x取值的间隔，这里先默认取0.5nm
length=0;
for j=1:layer
    length=length+thickness(j);
end
x=zeros(1,length/interval+1);
for j=1:length/interval+1
    x(j)=interval*(j-1);
end   %最后画图用得上x序列
nod=zeors(layer+1,1);  %创建一个节点数组，表示每个界面相应的x的“序列号”（而非相应x序列的值）
nod(1)=1;
for j=2:layer+1
    nod(j)=nod(j-1)+thickness(j-1)/interval;
end   %节点数组创建完毕
E=zeros(size(glass,1),length/interval+1);  %每一个波长对应一组相应x值的光电场分布
Q=zeros(size(glass,1),length/interval+1);  %每一个波长对应一组相应x值的能量吸收效率Q的分布
for lamda=1:size(glass,1)   %遍历所有波长，每个波长做出一个E(lamda，x)图和Q(lamda，x)图
    %我们要做准备工作，即用N把波长为glass（lamda，1）的相应所有复折射率存起来
    N=zeros(layer,1);
    for j=1:layer
        N(j)=epsilon{j}(lamda,2)+1i*epsilon{j}(lamda,3);
    end
    %进一步的准备工作是把每一层的ζ函数值存起来
    for j=1:layer
        Zetaelec(lamda,j)=2*pi*N(j)/glass(lamda,1);  %ζ函数公式是ζ=2πNj/λ
    end
    %首先，我们求出n层器件定义的n+1个界面传输矩阵Ijk
    %首先我们先求出玻璃基底和第一个薄层的界面I矩阵，以及最后一个薄层和出射空气层的界面I矩阵
    r=zeros(layer+1,1);
    t=zeors(layer+1,1);  %创建空的菲涅耳反射系数组和透射系数组
    Nglass=glass(lamda,2)+1i*glass(lamda,3);
    Nair=1;
    r(1)=(Nglass-N(1))/(Nglass+N(1));
    t(1)=2*N(1)/(Nglass+N(1));
    I{lamda,1}=(1/t(1))*[1,r(1);r(1),1];
    %这里求出了I(lamda,1)，也就是玻璃基底和第一个薄层的界面I矩阵
    r(layer+1)=(N(layer)-Nair)/(N(layer)+Nair);
    t(layer+1)=2*N(layer)/(N(layer)+Nair);
    I{lamda,layer+1}=(1/t(layer+1))*[1,r(layer+1);r(layer+1),1];
    %这里求出了I(lamda,layer+1)，也就是最后一个薄层和出射空气层的界面I矩阵
    for j=2:layer
        r(j)=(N(j-1)-N(j))/(N(j-1)+N(j));
        t(j)=2*N(j-1)/(N(j-1)+N(j));
        I{lamda,j}=(1/t(j))*[1,r(j);r(j),1];
    end
    %到这里为止，就求出了所有的层界面的I矩阵
    %接下来的任务是求取L矩阵（层矩阵）
    for j=1:layer
        L{lamda,j}=[exp(-1i*Zetaelec(lamda,j)*thickness(j)),0;0,exp(1i*Zetaelec(lamda,j)*thickness(j))];
    end
    %至此，L矩阵也求取完毕，接下来求取“任意点”的光电场值E（lamda，x）
    for j=1:layer  %首先遍历所有“层”
        left=nod(j);  %这是第j层的左边界对应的x序列的序列号，下同
        right=nod(j+1);
        if j<2           %求出第j层的左矩阵
            Sleft=I{lamda,1};
        else
            Sleft=I{lamda,1};
            for k=1:j-1
                Sleft=Sleft*L{lamda,k}*I{lamda,k+1};
            end
        end
        if j>layer-1     %求出第j层的右矩阵
            Sright=I{lamda,layer+1};
        else
            Sright=1;
            for k=j:layer-1
                Sright=Sright*I{lamda,k+1}*L{lamda,k+1};
            end
            Sright=Sright*I{lamda,layer+1};
        end
        Extinct=epsilon{j}(lamda,3); %第j层材料的消光系数
        nj=epsilon{j}(lamda,2);  %第j层材料的折射率
        for k=left:right   %求第j层中各个位置的光电场大小，以及相应的Q的大小
            tmp1=exp(-1i*Zetaelec(lamda,j)*thickness(j));
            tmp2=exp(1i*Zetaelec(lamda,j)*thickness(j));
            tmp3=exp(-1i*Zetaelec(lamda,j)*(thickness(j)-(k-left)*interval));
            tmp4=exp(1i*Zetaelec(lamda,j)*(thickness(j)-(k-left)*interval));    
            numerator=Sright(1,1)*tmp3+Sright(2,1)*tmp4;
            denominator=Sleft(1,1)*Sright(1,1)*tmp1+Sleft(1,2)*Sright(2,1)*tmp2;
            E(lamda,k)=E0in*numerator/denominator;
            Q(lamda,k)=2*pi*c*epsilon0*Extinct*nj*power(E(lamda,k),2)/glass(lamda,1);
        end
        subplot(1,2,1);
        plot(x,E(lamda,:));
        hold on;
        subplot(1,2,2);
        plot(x,Q(lamda,:));
        hold on;
    end
end

