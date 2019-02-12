c=299792458; %����
epsilon0=8.854187817*power(10,-12); %��ս�糣����0 = 8.854187817*10-12F/m
E0in=1;   %�������䵽��������һ�������Ĺ�糡�Ĵ�С��Ĭ����Ϊһ
layers=4;  %�����������ж��ٸ�����
epsilon=cell(layer,1);  %�������ڴ�����������纯�����������ʣ���Ԫ������
glass=load('glass.txt');  %��ȡ�����Ľ�纯��
for j=1:layer
    tmp=[num2str(j),'_epsilon.txt']; %��ȡ����Ľ�纯�����������ʣ�
    epsilon{j}=load(tmp);
end
thickness=load('thickness.txt');  %�������ĺ��ֵ
inmode=1; %����ģʽ��1��ʾ��ֱ���䣬2��ʾƫ��pƫ��3��ʾƫ��sƫ��
I=cell(size(glass,1),layer+1);
L=cell(size(glass,1),layer); 
Zetaelec=zeros(size(glass,1),layer);
%�ֱ𴴽��յĽ��������Ͳ�����飬�Լ������zeta����ֵ��������
%Ȼ������������Ĺ����䷽��������Ữ�ֶַ�ֵ����ȡ��x����
interval=0.5;  %�趨xȡֵ�ļ����������Ĭ��ȡ0.5nm
length=0;
for j=1:layer
    length=length+thickness(j);
end
x=zeros(1,length/interval+1);
for j=1:length/interval+1
    x(j)=interval*(j-1);
end   %���ͼ�õ���x����
nod=zeors(layer+1,1);  %����һ���ڵ����飬��ʾÿ��������Ӧ��x�ġ����кš���������Ӧx���е�ֵ��
nod(1)=1;
for j=2:layer+1
    nod(j)=nod(j-1)+thickness(j-1)/interval;
end   %�ڵ����鴴�����
E=zeros(size(glass,1),length/interval+1);  %ÿһ��������Ӧһ����Ӧxֵ�Ĺ�糡�ֲ�
Q=zeros(size(glass,1),length/interval+1);  %ÿһ��������Ӧһ����Ӧxֵ����������Ч��Q�ķֲ�
for lamda=1:size(glass,1)   %�������в�����ÿ����������һ��E(lamda��x)ͼ��Q(lamda��x)ͼ
    %����Ҫ��׼������������N�Ѳ���Ϊglass��lamda��1������Ӧ���и������ʴ�����
    N=zeros(layer,1);
    for j=1:layer
        N(j)=epsilon{j}(lamda,2)+1i*epsilon{j}(lamda,3);
    end
    %��һ����׼�������ǰ�ÿһ��Ħƺ���ֵ������
    for j=1:layer
        Zetaelec(lamda,j)=2*pi*N(j)/glass(lamda,1);  %�ƺ�����ʽ�Ǧ�=2��Nj/��
    end
    %���ȣ��������n�����������n+1�����洫�����Ijk
    %��������������������׺͵�һ������Ľ���I�����Լ����һ������ͳ��������Ľ���I����
    r=zeros(layer+1,1);
    t=zeors(layer+1,1);  %�����յķ���������ϵ�����͸��ϵ����
    Nglass=glass(lamda,2)+1i*glass(lamda,3);
    Nair=1;
    r(1)=(Nglass-N(1))/(Nglass+N(1));
    t(1)=2*N(1)/(Nglass+N(1));
    I{lamda,1}=(1/t(1))*[1,r(1);r(1),1];
    %���������I(lamda,1)��Ҳ���ǲ������׺͵�һ������Ľ���I����
    r(layer+1)=(N(layer)-Nair)/(N(layer)+Nair);
    t(layer+1)=2*N(layer)/(N(layer)+Nair);
    I{lamda,layer+1}=(1/t(layer+1))*[1,r(layer+1);r(layer+1),1];
    %���������I(lamda,layer+1)��Ҳ�������һ������ͳ��������Ľ���I����
    for j=2:layer
        r(j)=(N(j-1)-N(j))/(N(j-1)+N(j));
        t(j)=2*N(j-1)/(N(j-1)+N(j));
        I{lamda,j}=(1/t(j))*[1,r(j);r(j),1];
    end
    %������Ϊֹ������������еĲ�����I����
    %����������������ȡL���󣨲����
    for j=1:layer
        L{lamda,j}=[exp(-1i*Zetaelec(lamda,j)*thickness(j)),0;0,exp(1i*Zetaelec(lamda,j)*thickness(j))];
    end
    %���ˣ�L����Ҳ��ȡ��ϣ���������ȡ������㡱�Ĺ�糡ֵE��lamda��x��
    for j=1:layer  %���ȱ������С��㡱
        left=nod(j);  %���ǵ�j�����߽��Ӧ��x���е����кţ���ͬ
        right=nod(j+1);
        if j<2           %�����j��������
            Sleft=I{lamda,1};
        else
            Sleft=I{lamda,1};
            for k=1:j-1
                Sleft=Sleft*L{lamda,k}*I{lamda,k+1};
            end
        end
        if j>layer-1     %�����j����Ҿ���
            Sright=I{lamda,layer+1};
        else
            Sright=1;
            for k=j:layer-1
                Sright=Sright*I{lamda,k+1}*L{lamda,k+1};
            end
            Sright=Sright*I{lamda,layer+1};
        end
        Extinct=epsilon{j}(lamda,3); %��j����ϵ�����ϵ��
        nj=epsilon{j}(lamda,2);  %��j����ϵ�������
        for k=left:right   %���j���и���λ�õĹ�糡��С���Լ���Ӧ��Q�Ĵ�С
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

