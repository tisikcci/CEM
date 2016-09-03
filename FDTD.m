%%   ****  基于时域有限差分法对对称带状线加载矩形波导TM模研究   ****  %%
%%   ***************** 姓名：肖利 ******************       %%
%%
%基于FDTD法对对称屏蔽带状线TM模的计算分析
clear all
clc

global N rho mue epsilon c b a h M
%N为单边网格数，rho为a/b的值,mue为真空中的磁导率，epsiloon为真空中的介电常数，
%c为光速，b为边长的一半，a为微带长度的一半，h为步长，微带的厚度不计,M为微带单
%边划分网格数
time = cputime;         %cpu运行时间消耗
N = 30;%input('输入网格数(限为偶数，且不小于4)=');
c = 3e8;
mue = pi*4e-7;
epsilon = 1/(pi*36e9);
b = 5e-3;
h = 2*b/N;
rho = 0.7;%input('输入a/b的值(0.3-0.7之间)=');
a = rho*b;
M = round((b-a)/h);
%%
%....................................................................................
%构造矩阵
%构造I矩阵
I = eye(N-1);
%构造除微带线外A（i,j）对应的系数矩阵D,主要是其它对角矩阵块
D(1:N-1,1:N-1) = 0;
for i = 1:N-2;
    D(i,i+1) = 1;
    D(i+1,i) = 1;
end
for i = 1:N-1
    D(i,i) = -4;
end

%构造微带矩阵stripline，即为微带的贡献
for i = 1:N-1
    stripline(i,i) = -4;
end
for i = 2:N-2
    stripline(i,i-1) = 1;
    stripline(i,i+1) = 1;
end
i = N/2-M-1;
stripline(i,i-1) = 1;
stripline(i,i+1) = 0;
i = N/2+M+1;
stripline(i,i-1) = 0;
stripline(i,i+1) = 1;
for i = (N/2-M):(N/2+M)
    stripline(i,i-1) = 0;
    stripline(i,i+1) = 0;
end
stripline(1,2) = 1;
stripline(N-1,N-2) = 1;
%构造微带线上下两侧边界的网格的系数矩阵up_down,即A(i,j+1),A(i,j-1)的网格系数
for i = 1:N-1
    if i>=(N/2-M) && i<=(N/2+M)
        up_down(i,i) = 0;
    else
        up_down(i,i) = 1;
    end
end
%构造K矩阵
K(1:(N-1)*(N-1),1:(N-1)*(N-1)) = 0;
for i = 2:N-2
    if i == N/2-1
        K(((i-1)*(N-1)+1):(i*(N-1)),((i-2)*(N-1)+1):((i-1)*(N-1))) = I;
        K(((i-1)*(N-1)+1):(i*(N-1)),((i-1)*(N-1)+1):(i*(N-1))) = D;
        K(((i-1)*(N-1)+1):(i*(N-1)),(i*(N-1)+1):((i+1)*(N-1))) = up_down;
    elseif i == N/2;
             K(((i-1)*(N-1)+1):(i*(N-1)),((i-2)*(N-1)+1):((i-1)*(N-1))) = up_down;
             K(((i-1)*(N-1)+1):(i*(N-1)),((i-1)*(N-1)+1):(i*(N-1))) = stripline;
             K(((i-1)*(N-1)+1):(i*(N-1)),(i*(N-1)+1):((i+1)*(N-1))) = up_down;
    elseif i == N/2+1
                K(((i-1)*(N-1)+1):(i*(N-1)),((i-2)*(N-1)+1):((i-1)*(N-1))) = up_down;
                K(((i-1)*(N-1)+1):(i*(N-1)),((i-1)*(N-1)+1):(i*(N-1))) = D;
                K(((i-1)*(N-1)+1):(i*(N-1)),(i*(N-1)+1):((i+1)*(N-1))) = I;
    else
             K(((i-1)*(N-1)+1):(i*(N-1)),((i-2)*(N-1)+1):((i-1)*(N-1))) = I;
             K(((i-1)*(N-1)+1):(i*(N-1)),((i-1)*(N-1)+1):(i*(N-1))) = D;
             K(((i-1)*(N-1)+1):(i*(N-1)),(i*(N-1)+1):((i+1)*(N-1))) = I;  
    end
end
K(1:N-1,1:N-1) = D;
K(1:N-1,N:2*(N-1)) = I;
K((N-2)*(N-1) + 1:(N-1)^2,(N-3)*(N-1) + 1:(N-2)*(N-1)) = I;
K((N-2)*(N-1) + 1:(N-1)^2,(N-2)*(N-1) + 1:(N-1)^2) = D;
%%
%.........................................................................................
%相关值计算
[phi,value] = eig(K);%求出特征值向量矩阵phi和特征值eigv
eig_freq = c*sqrt(-value)/(2*pi*h);%根据特征值求出特征频率
fc = eig_freq;
fc(fc == 0) = inf;
fc_min1 = min(min(fc));
[x1,y1] = find(fc==min(min(fc)));
fc(:,y1) = inf;
fc_min2 = min(min(fc));
[x2,y2] = find(fc==min(min(fc)));
phi1 = phi(:,y1);
phi2 = phi(:,y2);
%.........................................................................................
%绘图
f = 30e9:0.25e9:40e9;
beta1 = real(2*pi*sqrt(mue*epsilon*(f.^2-fc_min1.^2)));
beta2 =real(2*pi*sqrt(mue*epsilon*(f.^2-fc_min2.^2)));
f = f/1e9;
plot(f,beta1,'rd--',f,beta2,'kp-')
axis normal
title('对称屏蔽带状线的相位常数(a/b=0.5,N=30)');
xlabel('频率（Hz）');ylabel('相位常数（rad/m）');
legend('最低TM模','次低TM模');
%..........................................................................................
 %%   
%场分布
kc1 = 2*pi*fc_min1/c;
kc2 = 2*pi*fc_min2/c;
w1 = 2*pi*fc_min1;
w2 = 2*pi*fc_min2;
Z0 = sqrt(mue/epsilon);
Ez1(1:N-1,1:N-1) = 0;
Ez2(1:N-1,1:N-1) = 0;
for n = 1:N-1
    Ez1(n,1:N-1) = phi1((n-1)*(N-1)+1:n*(N-1));
    Ez2(n,1:N-1) = phi2((n-1)*(N-1)+1:n*(N-1));
end
j = j;
for k = 1:N-1
    for i = 2:N-2
        Hy1(i,k) = (-j*w1*epsilon/(2*h*kc1^2))*(Ez1(i+1,k)-Ez1(i-1,k));
        Hy2(i,k) = (-j*w2*epsilon/(2*h*kc2^2))*(Ez2(i+1,k)-Ez2(i-1,k));
    end
    Hy1(1,k) = (-j*w1*epsilon/(2*h*kc1^2))*(Ez1(2,k));
    Hy2(1,k) = (-j*w2*epsilon/(2*h*kc2^2))*(Ez2(2,k));
    Hy1(N-1,k) = (j*w1*epsilon/(2*h*kc1^2))*(Ez1(N-2,k));
    Hy2(N-1,k) = (j*w2*epsilon/(2*h*kc2^2))*(Ez2(N-2,k));
end
for i = 1:N-1
    for k = 2:N-2
        Hx1(i,k) = (j*w1*epsilon/(2*h*kc1^2))*(Ez1(i,k+1)-Ez1(i,k-1));
        Hx2(i,k) = (j*w2*epsilon/(2*h*kc2^2))*(Ez2(i,k+1)-Ez2(i,k-1));
    end
    Hx1(i,1) = (j*w1*epsilon/(2*h*kc1^2))*(Ez1(i,2));
    Hx2(i,1) = (j*w2*epsilon/(2*h*kc2^2))*(Ez2(i,2));
    Hx1(i,N-1) = -(j*w1*epsilon/(2*h*kc1^2))*(Ez1(i,N-2));
    Hx2(i,N-1) = -(j*w2*epsilon/(2*h*kc2^2))*(Ez2(i,N-2));
end
Ex1 = Z0 * Hy1;
Ex2 = Z0 * Hy2;
Ey1 = -Z0 * Hx1;
Ey2 =  -Z0 * Hx2;
axis([0 40 0 40])
subplot(2,2,1);
Efield1 = quiver(imag(Ey1),imag(Ex1));%第一个TM模的横截面电场分布图
title('第一个TM模的横截面电场分布图')
subplot(2,2,2);
Hfield1 = quiver(imag(Hy1),imag(Hx1),'r');%第一个TM模的横截面磁场分布图
title('第一个TM模的横截面磁场分布图')
subplot(2,2,3);
Efield2 = quiver(imag(Ey2),imag(Ex2));%第二个TM模的横截面电场分布图
axis([0 30 0 30])
title('第二个TM模的横截面电场分布图')
subplot(2,2,4);
Hfield2 = quiver(imag(Hy2),imag(Hx2),'r');%第二个TM模的横截面磁场分布图
title('第二个TM模的横截面磁场分布图') 

