%书上189页算例
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


for t=1:1:m+1
    a(t,1)=0;
    a(1,t)=0;
    a(t,m+1)=0;
    a(m+1,t)=0;
end

a(1,1)=1;
a(m+1,m+1)=1;
b(1,1)=0;
b(m+1,1)=0;

u=a\b;

function [z] = f(x)
    z=((16*pi^2-4)*sin(4*pi*x)-16*pi*cos(4*pi*x))*exp(2*x);
end


function [z] = phi(i,x)
global xco h;
temp=abs(x-xco(i,1));
if(temp<=h)
    z=1.0-temp/h;
else
    z=0;
end
end


function [z] = fun1(i,x)
    z=f(x)*phi(i,x);
end


function [z] = intergral(a,b,i)
    gauss=0.5773502692;
    mid=(a+b)/2;
    w=(b-a)/2;
    z=w*(fun1(i,mid+w*gauss)+fun1(i,mid-w*gauss));
end



