%205页案例,抛物型方程

clc;clear
tic
global xco h;    %全局变量
m=32;             %空间剖分数
n=6400;             %时间剖分数
elem=m;          %单元数
h=pi/m;          %   空间步长
T=1.0;           %   时间终点
tau=T/n;         %   时间步长
a1=2.0;          %   常系数a
alpha=zeros(m+1,n+1); %存储所有解
CN=0.5;      %CN法时间推进


xco=zeros(m+1,1);   %存储x方向节点的位置，共计m+1个节点
for i=1:1:m+1
    xco(i,1)=(i-1)*h;
end

t=zeros(n+1,1);   %存储t方向节点的位置，共计n+1个节点
for i=1:1:n+1
    t(i,1)=(i-1)*tau;
end

lnd=zeros(elem,2);  %第一个索引是单元，第二个索引是局部节点编号，存储第几个单元的第几个节点的全局编号是多少
for i=1:1:m
    lnd(i,1)=i;
    lnd(i,2)=i+1;
end

% 系数矩阵A、B不变，先组装好
A=zeros(m+1,m+1);   %A是phii与phij的积分
B=zeros(m+1,m+1);   %A是phii导与phij导的积分

Aelem1=h/3;
Aelem2=h/6;
Aelem=[Aelem1,Aelem2;Aelem2,Aelem1];    %A单元系数矩阵

Belem_1=zeros(2,1);
Belem_1(1,1)=-1.0/h;   %形函数导数是正负1/h，
Belem_1(2,1)=1.0/h;

Belem=zeros(2,2);     %B每个单元上系数矩阵，在一个单元上相邻或同一形函数积分那就是(1/h)^2*h=1/h
for i=1:1:2
for j=1:1:2
    Belem(i,j)=Belem_1(i,1)*Belem_1(j,1)*h;
end
end

for e=1:1:elem %遍历单元   
    for i=1:1:2
    for j=1:1:2
        row=lnd(e,i);
        coln=lnd(e,j);
        A(row,coln)=A(row,coln)+Aelem(i,j);  %组装整体系数矩阵
        B(row,coln)=B(row,coln)+Belem(i,j);  %组装整体系数矩阵
    end
    end
end

%初始解
for i=1:1:m+1
    alpha(i,1)=psi(xco(i,1));
end  

d=zeros(m+1,1); %右端项F 
AA=zeros(m+1,m+1);%最终计算的系数矩阵AA
AA=A+CN*a1*tau*B;
BB1=zeros(m+1,1);


% 修改边界条件
AA(1,1)=1;
AA(2,1)=0;
AA(1,2)=0;
AA(m+1,m+1)=1;
AA(m+1,m)=0;
AA(m,m+1)=0;

for k=1:1:n     %时间循环
    %组装右端源项f
    tmid=(t(k)+t(k+1))/2;
    for i=2:m
        d(i,1)=intergral(xco(i-1,1),xco(i,1),i,tmid)+intergral(xco(i,1),xco(i+1,1),i,tmid);
    end
    d(1,1)=intergral(xco(1,1),xco(2,1),1,tmid);
    d(m+1,1)=intergral(xco(m,1),xco(m+1,1),m+1,tmid);

    BB1=(A-(1-CN)*a1*tau*B)*alpha(:,k);
    BB=BB1+tau*d;

%     修改边界条件
    BB(1,1)=0;
    BB(m+1,1)=0;
    alpha(:,k+1)=AA\BB;
end

result=zeros(m+1,n+1);
for i=1:1:m+1
    for j=1:1:n+1
        result(i,j)=exact(xco(i,1),t(j,1));
    end
end

sum(sum(result-alpha))/n
% for i=1:1:n+1
%     plot(result(:,i));
%     hold on;
%     plot(alpha(:,i));
%     legend('真解','初始误差较小的有限元解')
%     pause(0.1)
%     clf;
% end




function [z] = f(x,t)       %方程的f,右端项，源项
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

function [z]=psi(x)   %初始条件
    z=5*sin(2*x)-30*sin(3*x)+2;
%           z=1;
%       end
end


function [z] = fun1(i,x,t)   %右端荷载
    z=f(x,t)*phi(i,x);
end


function [z] = intergral(a,b,i,t)   %在区间ab上对右端荷载采用两点高斯公式积分
    gauss=0.5773502692;
    mid=(a+b)/2;
    w=(b-a)/2;
    z=w*(fun1(i,mid+w*gauss,t)+fun1(i,mid-w*gauss,t));
end

function [z]=exact(x,t)
    z=5*exp(-8.0*t)*sin(2*x)-30*exp(-18.0*t)*sin(3*x);
end


