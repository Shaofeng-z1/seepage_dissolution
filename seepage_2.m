%可设置第一类边界条件

clc;clear
global xco h;
m=8;
elem=m;
h=1.0/m;


xco=zeros(m+1,1);
for i=1:1:m+1
xco(i,1)=(i-1)*h;
end

lnd=zeros(elem,2);
for i=1:1:m
    lnd(i,1)=i;
    lnd(i,2)=i+1;
end

a=zeros(m+1,m+1);
b=zeros(m+1,1);

alpha=zeros(2,1);
alpha(1,1)=-1.0/h;
alpha(2,1)=1.0/h;

ea=zeros(2,2);
for i=1:1:2
for j=1:1:2
ea(i,j)=alpha(i,1)*alpha(j,1)*h;
end
end

g=zeros(2,1);
for e=1:1:elem
    i=lnd(e,1);
    j=lnd(e,2);
    for k=1:1:2
        g(k,1)=intergral(xco(i,1),xco(j,1),lnd(e,k));
    end
    
    for i=1:1:2
    for j=1:1:2
        row=lnd(e,i);
        coln=lnd(e,j);
        a(row,coln)=a(row,coln)+ea(i,j);
    end
        k=lnd(e,i);
        b(k)=g(i)+b(k);
    end
end


%修改边界条件
big=10^8;
a(1,1)=big;
a(8,8)=big;
a(m+1,m+1)=big;
b(1,1)=a(1,1)*10;
b(8,1)=a(8,8)*20;
b(m+1,1)=a(m+1,m+1)*0;

u=a\b;

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



