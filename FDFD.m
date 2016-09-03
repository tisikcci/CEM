%%   ****  基于频域有限差分法对对称带状线加载矩形波导TM模研究   ****  %%
%%   ***************** 姓名：肖利 ******************       %%
%%   *************  学号：201421040239  *************     %%
%%模型示意图如下
%         *****************************       ——
%         *                           *
%         *                           *        宽
%         *     |    导带宽w   |      *        by
%         *     ************** |缝宽a|*  ——  
%         *                           *   高     
%         *                           *   h
%         *****************************  ——  ——                      
%         |          波导长bx          |
%% 基本数据处理
clc;clear;close all;
b=input('请输入波导宽度一半b的值(单位：mm)：'); 
b=b*1e-3; % b=5mm  宽的一半
by=2*b;  %%波导的宽 
bx=by;  %%波导的长,算例1情况。
rou=input('请输入a/b的值(0.3到0.7之间)：'); 
a=b*rou;%由rou=a/b,计算a值
w=bx-2*a;%带状线的高度及宽度
h=1/2*by;    %导带高          
epsilon = 8.85e-12; %真空中的介电常数
mu=4*pi*(1e-7); %真空中的磁导率
epsilonxyz=1;%相对介电常数分量
c=3e8;%真空中光速
jim=sqrt(-1);%虚数i用jim（j-image)表示
%% 划分网格x方向nx个，y方向ny个
nx=input('请输入x方向划分的网格数：');nx=ceil(nx);
ny=input('请输入y方向划分的网格数：');ny=ceil(ny);disp(sprintf('划分网格数：x方向 %g 个，y方向 %g 个',nx,ny))
tic;  %%开始计时
ni=nx+1;nj=ny+1;%x与y方向结点数目
dx=bx/nx;dy=by/ny;%x方向与y方向网格步长
s1=ceil(a/dx)+1;
s2=floor((a+w)/dx)+1;
sy=round(h/dy)+1;   %计算导带位置节点编号
%% 提取系数矩阵
%%四个系矩阵均为 toeplitz 矩阵
B1(ni*ny,ni*nj)=0;%%EZ 产生HX。由HX,EZ控制行列大小
B2(nx*nj,ni*nj)=0;%%EZ 产生Hy。由HY,EZ控制行列大小
B3(ni*nj,ni*ny)=0;%%HX 产生部分EZ。由EZ,HX控制行列大小
B4(ni*nj,nx*nj)=0;%%HY 产生部分EZ。由EZ,HY控制行列大小
B(ni*ny+nx*nj+ni*nj,ni*ny+nx*nj+ni*nj)=0;%%B1+B2+B3行和，B3+B4+B1列和
%% hx ni*ny个 hy nx*nj 个 EZ ni*nj个
c1(ny)=0;c1(1)=-1;        %构造TOEPLITZ矩阵B1_sub所需列
r1(nj)=0;r1(1)=-1;r1(2)=1; %构造TOEPLITZ矩阵B1_sub所需行
B1_sub=toeplitz(c1,r1);   %B1
for i=1:ni
    B1((i-1)*ny+(1:ny),(i-1)*nj+(1:nj))=B1_sub(1:ny,1:nj);
end
B1=jim/dy*B1;
c2(nx*nj)=0;c2(1)=-1;        %构造TOEPLITZ矩阵B2所需列
r2(ni*nj)=0;r2(1)=-1;r2(1+nj)=1; %构造TOEPLITZ矩阵B2所需行
B2=-jim/dx*toeplitz(c2,r2);   %B2
c3(nj)=0;c3(1)=1;c3(2)=-1;        %构造TOEPLITZ矩阵B3_sub所需列
r3(ny)=0;r3(1)=1; %构造TOEPLITZ矩阵B3_sub所需行
B3_sub=toeplitz(c3,r3); 
for i=1:ni
    B3((i-1)*nj+(1:nj),(i-1)*ny+(1:ny))=B3_sub(1:nj,1:ny);
end
B3=jim/dy/epsilonxyz*B3;          %B3
c4(ni*nj)=0; c4(1)=1; c4(1+nj)=-1 ;     %构造TOEPLITZ矩阵B4所需列
r4(nx*nj)=0;r4(1)=1;  %构造TOEPLITZ矩阵B4所需行
B4=-jim/dx/epsilonxyz*toeplitz(c4,r4);   %B4
B(1:ni*ny,ni*ny+nx*nj+(1:ni*nj))=B1(1:ni*ny,1:ni*nj);
B(ni*ny+(1:nx*nj),ni*ny+nx*nj+(1:ni*nj))=B2(1:nx*nj,1:ni*nj);
B(ni*ny+nx*nj+(1:ni*nj),1:ni*ny)=B3(1:ni*nj,1:ni*ny);
B(ni*ny+nx*nj+(1:ni*nj),ni*ny+(1:nx*nj))=B4(1:ni*nj,1:nx*nj);
%边界处理 Ez=0
B(ni*ny+nx*nj+(1:nj),:)=0; %%左边界
B(ni*ny+nx*nj+ni*nj-nj+(1:nj),:)=0; %%右边界
for i=1:ni
B(ni*ny+nx*nj+1+nj*(i-1),:)=0 ;  %%下边界
B(ni*ny+nx*nj+i*nj,:)=0 ;   %%上边界
end
%%加入金属导带
for i=s1:s2
    B(ni*ny+nx*nj+sy+nj*(s1-1:s2-1),:)=0;
end

%% 计算系数矩阵特征值与特征向量
%mode=input('输入要得到的TM模式截止频率数：');
mode=5;
[EH, k0]=eig(B);
%% 一、计算截止频率
k0=real(k0);
f=k0*c/(2*pi);
f_sort=sort(f(:));
fc(1,1:mode)=0;
ii = 0;
%删除伪解
for i = 1:length(f_sort)
if f_sort(i) > 1e9
ii = ii+1;
fc(ii) = f_sort(i);
if mode== ii
break;
end
end
end
disp(['运行时间：',num2str(toc)]);
%输出截止频率
disp('***********************************************');
disp('FDFD计算结果：');
for i=1:mode
disp(sprintf('TM波第 %g 次模截止频率：fc= %gGHz；',i,fc(i)/(1e9)));
end
%% 二、计算场分量
[x1, y1] = find(f == fc(1));%得到第一个截止频率所对应的场分量 f是矩阵X1,Y1表FC(1)所在行，列。
[x2, y2] = find(f == fc(2));%得到第二个截止频率所对应的场分量
% 第一个模式的场分量
HX1=EH(1:ni*ny, y1);%提取Hx分量
HX1=reshape(HX1, ny, ni);%转化为矩阵形式
if  sum(sum(abs(real(HX1)))) > sum(sum(abs(imag(HX1))))  
    HX1 = real(HX1);
else
    HX1 = imag(HX1);
end
HY1 = EH(ni*ny+(1:nx*nj), y1);%提取Hy分量
HY1 = reshape(HY1,nj, nx);%转化为矩阵形式
if  sum(sum(abs(real(HY1)))) > sum(sum(abs(imag(HY1))))  
    HY1 = real(HY1);
else
    HY1 = imag(HY1);
end
EZ1 = EH(ni*ny+nx*nj+(1:ni*nj), y1);%提取Ez分量
EZ1 = reshape(EZ1,nj, ni);%转化为矩阵形式
if  sum(sum(abs(real(EZ1)))) > sum(sum(abs(imag(EZ1))))  
    EZ1 = real(EZ1);
else
    EZ1 = imag(EZ1);
end
% 第二个模式的场分量
HX2=EH(1:ni*ny, y2);%提取Hx分量
HX2=reshape(HX2, ny, ni);%转化为矩阵形式
if  sum(sum(abs(real(HX2)))) > sum(sum(abs(imag(HX2))))  
    HX2 = real(HX2);
else
    HX2 = imag(HX2);
end
HY2 = EH(ni*ny+(1:nx*nj), y2);%提取Hy分量
HY2 = reshape(HY2,nj, nx);%转化为矩阵形式
if  sum(sum(abs(real(HY2)))) > sum(sum(abs(imag(HY2))))  
    HY2 = real(HY2);
else
    HY2 = imag(HY2);
end
EZ2 = EH(ni*ny+nx*nj+(1:ni*nj), y2);%提取Ez分量
EZ2 = reshape(EZ2,nj, ni);%转化为矩阵形式
if  sum(sum(abs(real(EZ2)))) > sum(sum(abs(imag(EZ2))))  
    EZ2 = real(EZ2);
else
    EZ2 = imag(EZ2);
end
%% 画场分量三维图
figure
subplot(1,3,1)
surf(HX1);
title('TM模1 Hx分量X-Y平面视图')
%figure
subplot(1,3,2)
surf(HY1);
title('TM模1 Hy分量X-Y平面视图')
subplot(1,3,3)
%figure
surf(EZ1);
title('TM模1 Ez分量X-Y平面视图')
figure
subplot(1,3,1)
surf(HX2);
title('TM模2 Hx分量X-Y平面视图')
%figure
subplot(1,3,2)
surf(HY2);
title('TM模2 Hy分量X-Y平面视图')
%figure
subplot(1,3,3)
surf(EZ2);
title('TM模2 Ez分量X-Y平面视图')
%% 绘制色散特性曲线
f1=fc(1):0.5e9:50e9;
f2=fc(2):0.5e9:50e9;
beta1 = 2*pi*sqrt(epsilon*mu)*sqrt(f1.^2 -(fc(1))^2);
beta2 = 2*pi*sqrt(epsilon*mu)*sqrt(f2.^2 -(fc(2))^2);
beta1=abs(beta1);beta2=abs(beta2);
f1=f1/(1e9);  %单位为GHz
f2=f2/(1e9);  %单位为GHz
figure
plot(f1,beta1,'-*r');  %第一高次模
hold on;
plot(f2,beta2,'-pb');  %第二高次模
hold off;
ylabel('rad/m');
xlabel('频率(GHz)');
ylabel('相位常数（rad/m）')
title('非对称屏蔽带状线矩形波导TM模的色散特性')
legend('TM第一高次模','TM第二高次模','Location','NorthWest')
%% 画横截面场分量二维分布图
%%根据横纵关系计算电场切向分量
f1_g =30e9;  %%设置模式1的一个工作频率30GHZ
beta1_g=2*pi*sqrt(epsilon*mu)*sqrt(f1_g^2 -(fc(1))^2);
omiga1=2*pi*f1_g;
k1=2*pi*f1_g/c;
ZTM1=omiga1*mu*k1/beta1_g^2;
HX_1(1:ny, 1:nx) = (HX1(1:ny, 1:nx) + HX1(1:ny, 2:ni))/2;
HY_1(1:ny, 1:nx) = (HY1(1:ny, 1:nx) + HY1(2:nj, 1:nx))/2;
EX_1= ZTM1*HY_1;
EY_1= -ZTM1*HX_1;
%figure
subplot(2,2,1)
quiver(EX_1, EY_1);
title('TM第1模式 横截面电场分量二维图')
subplot(2,2,2)
quiver(HX_1, HY_1,'r');
title('TM第1模式 横截面磁场分量二维图')
f2_g = 40e9;  %%设置模式2的一个工作频率40GHZ
beta2_g=2*pi*sqrt(epsilon*mu)*sqrt(f2_g^2 -(fc(2))^2);
omiga2=2*pi*f2_g;
k2=2*pi*f2_g/c;
ZTM2=omiga2*mu*k2/beta2_g^2;
HX_2(1:ny, 1:nx) = (HX2(1:ny, 1:nx) + HX2(1:ny, 2:ni))/2;
HY_2(1:ny, 1:nx) = (HY2(1:ny, 1:nx) + HY2(2:nj, 1:nx))/2;
EX_2= ZTM2*HY_2;
EY_2= -ZTM2*HX_2;
%Figure 
subplot(2,2,3)
quiver(EX_2, EY_2);
title('TM第2模式 横截面电场分量二维图')
subplot(2,2,4)
quiver(HX_2, HY_2,'r');
title('TM第2模式 横截面磁场分量二维图')

