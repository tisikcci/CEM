%---------------------------PO_PEC双站实现-------------------------%
%2016-4-12
%
clear;
clc;
p = load('NODE.txt');
t = load('FACE.txt');
p(:,1)=[];   %删除第一列
t(:,1)=[];   %删除第一列
TrianglesTotal=length(t);          %总的三角形面元个数
%%
%每个面元的面积、重心、外法矢
A = p(t(:,1),:);                   %每个三角面元的节点坐标
B = p(t(:,2),:);
C = p(t(:,3),:);
AB = B-A;                          %三角面元边矢
BC = C-B;
CA = A-C;
ABxBC = cross(AB,BC);
Area = 0.5*sqrt(sum(ABxBC.^2,2));  %每个面元的面积
Center = 1/3*(A+B+C);              %第m个三角面元的几何中心坐标(图心)
n_i = ABxBC./repmat(sqrt(sum(ABxBC.^2,2)),1,3);       %面元法矢
ToT_S_s = sum(Area);                                  %所有三角面元求和
r_n = [AB,BC,CA];                  %一行九列的三边矢量矩阵      
r_c = [(A+B)/2,(B+C)/2,(A+C)/2];   %一行九列的三边中点坐标矩阵
clear ABxBC;
%%
%电磁参数设置
f = 3e8;                %入射波频率设置  
epsilon_ = 8.854e-012;  %自由空间介电常数
mu_ = 1.257e-006;       %自由空间磁导率
c_=1/sqrt(epsilon_*mu_);%自由空间光速
lambda =c_/f;           %入射波波长
Z0=sqrt(mu_/epsilon_);  %本征阻抗
k0 = 2*pi/lambda;       %空间自由波数
r = 2e2*lambda;       %空间场点距离（自行设置，但要足够远）    
%%
%外加场设置
thta = 35;
de_thta = thta*pi/180;
phi = 30;
de_phi = phi*pi/180;
thta_s = 0:360;
AA = length(thta_s);
de_thta_s = thta_s'*pi/180;
phi_s = 30;
de_phi_s = phi_s*pi/180;
%设置电磁场的极化方向和波数方向
k_i= -[sin(de_thta)*cos(de_phi) sin(de_thta)*sin(de_phi) cos(de_thta)];          %入射波方向
e_i_theta = [cos(de_phi)*cos(de_thta) cos(de_thta)*sin(de_phi) -sin(de_thta)];  %thta极化入射
e_i_phi = [-sin(de_phi) cos(de_phi) 0];                             %fire极化入射
pol_e_i = e_i_theta;
h_i = cross(k_i,pol_e_i);   %入射磁场的极化方向
k_s = [sin(de_thta_s)*cos(de_phi_s) sin(de_thta_s)*sin(de_phi_s) cos(de_thta_s)];%散射波方向矢量
e_r_thta_s = [cos(de_thta_s)*cos(de_phi_s) cos(de_thta_s)*sin(de_phi_s) -sin(de_thta_s)];%thta极化接收
e_r_fire_s = repmat([-sin(de_phi_s) cos(de_phi_s) 0],AA,1);  %fire极化接收
%入射场设置
H0 = 1/Z0;     %入射场选为磁场，幅度选为1
%%
%做面元遮蔽判断：
%自遮蔽判断
afa0 = n_i*k_i';              %所有入射方向矢量与面元法矢作用(每一列即是对应一个角度下时，所有面片的遮蔽判断角)
sita = afa0<0;                %afa0<0的所有面元都置为1，afa0>=0的面元都被自身遮蔽，置为0
%互遮蔽判断
%此处预留再处理
clear afa0
%%
%计算任意远场空间处的散射场
%参数准备
part1 = 1j*k0*Z0*H0*exp(-1j*k0*r)/(2*pi*r);
part3 = -1j/k0;
W = repmat(k_i,AA,1)-k_s;
E_n=zeros(TrianglesTotal,3);
E = zeros(AA,3);
deta = lambda*1e-5;
for i=1:AA
    k_s_n = repmat(k_s(i,:),TrianglesTotal,1);
    h_i_n = repmat(h_i,TrianglesTotal,1);
    W_n = repmat(W(i,:),TrianglesTotal,1);
    W_i = W(i,:)';
    cross_W_n=cross(W_n,n_i);
    part2 = repmat(sita,1,3).*cross(k_s_n,cross(k_s_n,cross(n_i,h_i_n)));
    you1 = part1*part2;
    part4_1=sum(r_n(:,1:3).*cross_W_n,2).*exp(-1j*k0*r_c(:,1:3)*W_i).*sinc(k0/2*r_n(:,1:3)*W_i);
    part4_2=sum(r_n(:,4:6).*cross_W_n,2).*exp(-1j*k0*r_c(:,4:6)*W_i).*sinc(k0/2*r_n(:,4:6)*W_i);
    part4_3=sum(r_n(:,7:9).*cross_W_n,2).*exp(-1j*k0*r_c(:,7:9)*W_i).*sinc(k0/2*r_n(:,7:9)*W_i);
    you2 = part4_1+part4_2+part4_3;
    you3 = Area.*exp(-1j*k0*A*W_i);
    T_2 = sum(cross_W_n.^2,2);
    for j=1:TrianglesTotal
        if T_2(j,1)<deta;
            E_n(j,:)=you1(j,:)*you3(j,1);
        else
            E_n(j,:)=you1(j,:)*you2(j,1)*part3/T_2(j,1);
        end
    end
    E(i,:)=sum(E_n);
end
%%
%根据散射场的值计算RCS
rcs_thta = 4*pi*r^2*(abs(sum(E.*e_r_thta_s,2).^2));
rcs_phi = 4*pi*r^2*(abs(sum(E.*e_r_fire_s,2).^2));
RCS_THTA = 10*log10(rcs_thta);
RCS_PHI = 10*log10(rcs_phi);
%%
%绘图
figure(1)
plot(thta_s,RCS_THTA)
hold on
figure(2)
plot(thta_s,RCS_PHI)