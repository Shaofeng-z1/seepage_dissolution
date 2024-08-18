%注释完整

clc;clear
global xco h; %全局变量
m=100;          %剖分数
elem=m;       %单元数
h=1.0/m;      %步长
b_l=10000;       %左边界
b_r=0;        %右边界
big=10^8;     %修改边界条件需要的一个大数

xco=zeros(m+1,1);   %存储节点的位置，共计m+1个节点
for i=1:1:m+1
xco(i,1)=(i-1)*h;
end

lnd=zeros(elem,2);  %第一个索引是单元，第二个索引是局部节点编号，存储第几个单元的第几个节点的全局编号是多少
for i=1:1:m
    lnd(i,1)=i;
    lnd(i,2)=i+1;
end

a=zeros(m+1,m+1);    %系数矩阵
b=zeros(m+1,1);      %系数矩阵的右端项

alpha=zeros(2,1);    
alpha(1,1)=-1.0/h;   %形函数导数是正负1/h，
alpha(2,1)=1.0/h;

ea=zeros(2,2);     %每个单元上系数矩阵，在一个单元上相邻或同一形函数积分那就是(1/h)^2*h=1/h
for i=1:1:2
for j=1:1:2
ea(i,j)=alpha(i,1)*alpha(j,1)*h;
end
end

g=zeros(2,1);
for e=1:1:elem %遍历单元
    i=lnd(e,1);
    j=lnd(e,2);
    for k=1:1:2   %遍历单元上每个节点
        g(k,1)=intergral(xco(i,1),xco(j,1),lnd(e,k));  %计算每个单元上的节点的右端项
    end
    
    for i=1:1:2
    for j=1:1:2
        row=lnd(e,i);
        coln=lnd(e,j);
        a(row,coln)=a(row,coln)+ea(i,j);  %组装整体系数矩阵
    end
        k=lnd(e,i);
        b(k,1)=g(i,1)+b(k,1);             %组装整体右端项
    end
end


%修改边界条件
a(1,1)=big;
a(m+1,m+1)=big;
b(1,1)=a(1,1)*b_l;
b(m+1,1)=a(m+1,m+1)*b_r;

u=a\b;
plot(u);
function [z] = f(x)       %方程的f
    z=0;
end


function [z] = phi(i,x)   %第i个单元的行函数
global xco h;
temp=abs(x-xco(i,1));
if(temp<=h)
    z=1.0-temp/h;
else
    z=0;
end
end


function [z] = fun1(i,x)   %右端荷载
    z=f(x)*phi(i,x);
end


function [z] = intergral(a,b,i)   %在区间ab上对右端荷载采用两点高斯公式积分
    gauss=0.5773502692;
    mid=(a+b)/2;
    w=(b-a)/2;
    z=w*(fun1(i,mid+w*gauss)+fun1(i,mid-w*gauss));
end



