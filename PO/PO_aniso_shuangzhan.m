%-----------------PO各项异性材料双站计算实现-------------------%
%2016-4-18
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
x = AB./repmat(sqrt(sum(AB.^2,2)),1,3);
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
r = 2e3*lambda;       %空间场点距离（自行设置，但要足够远）  
%%
%目标体参数设置
afa = 0*pi/180;                   %各向异性材料主方向角偏角
Zuu = 0;%(0.1-0.01j)*Z0;%(0.1133+0.0166j)*Z0;        %主方向值设置
Zvv = 0;%(0.1-0.01j)*Z0;%(0.1133+0.0166j)*Z0;        %主方向值设置
Zxx = (cos(afa))^2*Zuu+(sin(afa))^2*Zvv;           %阻抗分量设置
Zxy = cos(afa)*sin(afa)*(Zuu-Zvv); 
Zyx = cos(afa)*sin(afa)*(Zuu-Zvv);
Zyy = (cos(afa))^2*Zvv+(sin(afa))^2*Zuu;
%%
%外加场设置
thta = 0;
de_thta = thta*pi/180;
phi = 0;
de_phi = phi*pi/180;
thta_s = 0:180;
AA = length(thta_s);
de_thta_s = thta_s'*pi/180;
phi_s = 0;
de_phi_s = phi_s*pi/180;
%设置电磁场的极化方向和波数方向
k_i= -[sin(de_thta)*cos(de_phi) sin(de_thta)*sin(de_phi) cos(de_thta)];          %入射波方向
e_i_theta = [cos(de_phi)*cos(de_thta) cos(de_thta)*sin(de_phi) -sin(de_thta)];  %thta极化入射
e_i_phi = [-sin(de_phi) cos(de_phi) 0];                             %fire极化入射
h_i_phi = cross(k_i,e_i_theta);   %入射磁场的垂直极化方向
h_i_theta = cross(k_i,e_i_phi);       %入射磁场的水平极化方向
k_s = [sin(de_thta_s)*cos(de_phi_s) sin(de_thta_s)*sin(de_phi_s) cos(de_thta_s)];%散射波方向矢量
e_r_thta_s = [cos(de_thta_s)*cos(de_phi_s) cos(de_thta_s)*sin(de_phi_s) -sin(de_thta_s)];%thta极化接收
e_r_fire_s = repmat([-sin(de_phi_s) cos(de_phi_s) 0],AA,1);  %fire极化接收
%入射场设置
H0 = 1;     %入射场选为磁场，幅度选为1
E0 = Z0*H0;  %电场幅值设置
%%
%根据入射场计算出每个面元的反射系数及场分量
%引入极化角概念：极化角为电场实际极化方向矢量与入射平面的夹角pol_j
pol_j = 0/180*pi;
E_thta = E0*cos(pol_j);   %与相位项无关的水平极化电场幅值表达
E_phi = E0*sin(pol_j);    %与相位项无关的垂直极化电场幅值表达
H_thta = H0*sin(pol_j);   %与相位项无关的水平极化磁场幅值表达
H_phi = H0*cos(pol_j);    %与相位项无关的垂直极化磁场幅值表达  
E_inc = E_thta*e_i_theta+E_thta*e_i_phi; %总坐标系电场矢量极化表达

%反射系数（与入射角和面片的位置有关即由：入射方向矢量和面片法向矢量确定）
e_ph = repmat(e_i_phi,TrianglesTotal,1);
e_tao = cross(n_i,cross(n_i,e_ph));
cos_i = -n_i*k_i';
sin_i = sqrt(1-cos_i.^2);
cos_X = sum(e_tao.*x,2);
sin_X = sqrt(1-cos_X.^2);
%用来求解反射系数的参数
Aa = 1/Z0*(Zxx*sin_X.^2+Zyy*cos_X.^2-Zxy*(sin_X.*cos_X)-Zyx*(sin_X.*cos_X));
Bb = 1/Z0*(Zxy*sin_X.^2+Zyx*cos_X.^2+Zxx*(sin_X.*cos_X)-Zyy*(sin_X.*cos_X));
Cc = 1/Z0*(Zxy*sin_X.^2-Zyx*cos_X.^2+Zxx*(sin_X.*cos_X)-Zyy*(sin_X.*cos_X));
Dd = -1/Z0*(Zyy*sin_X.^2+Zxx*cos_X.^2+Zxy*(sin_X.*cos_X)+Zyx*(sin_X.*cos_X));
%反射系数
R11 = (Dd.*cos_i.^2-Aa.*Dd.*cos_i+Bb.*Cc.*cos_i-cos_i+Aa)./(Dd.*cos_i.^2+Aa.*Dd.*cos_i-Bb.*Cc.*cos_i-cos_i-Aa);
R12 = (-2*Bb.*cos_i)./(Dd.*cos_i.^2+Aa.*Dd.*cos_i-Bb.*Cc.*cos_i-cos_i-Aa);
R21 = (-2*Cc.*cos_i)./(Dd.*cos_i.^2+Aa.*Dd.*cos_i-Bb.*Cc.*cos_i-cos_i-Aa);
R22 = (Dd.*cos_i.^2+Aa.*Dd.*cos_i-Bb.*Cc.*cos_i+cos_i+Aa)./(Dd.*cos_i.^2+Aa.*Dd.*cos_i-Bb.*Cc.*cos_i-cos_i-Aa);
%%
%做面元遮蔽判断：
%自遮蔽判断
afa0 = n_i*k_i';              %所有入射方向矢量与面元法矢作用(每一列即是对应一个角度下时，所有面片的遮蔽判断角)
sita = afa0<0;                %afa0<0的所有面元都置为1，afa0>=0的面元都被自身遮蔽，置为0
%互遮蔽判断
%此处预留再处理
clear afa0
%%
%散射体表面感应电流和感应磁流
J_part1 = ((1-R11).*E_phi-R21.*E_thta).*cos_i;
J_part2 = (R12.*E_phi+(1+R11).*E_thta);
Js = 1/Z0*(repmat(sita,1,3).*(repmat(J_part1,1,3).*e_ph+repmat(J_part2,1,3).*cross(n_i,e_ph)));
Jm_part1 = ((1-R11).*E_thta-R12.*E_phi).*cos_i;
Jm_part2 = R21.*E_thta+(1+R22).*E_phi;
Jm = repmat(sita,1,3).*(repmat(Jm_part1,1,3).*e_ph-repmat(Jm_part2,1,3).*cross(n_i,e_ph));
%%
%计算任意远场空间处的散射场
%参数准备
part1 = 1j*k0*Z0*exp(-1j*k0*r)/(4*pi*r);
part3 = -1j/k0;
W = repmat(k_i,AA,1)-k_s;
E_n=zeros(TrianglesTotal,3);
E = zeros(AA,3);
deta = lambda*1e-5;
for i=1:AA
    k_s_n = repmat(k_s(i,:),TrianglesTotal,1);
    W_n = repmat(W(i,:),TrianglesTotal,1);
    W_i = W(i,:)';
    cross_W_n=cross(W_n,n_i);
    part2_1 = cross(k_s_n,cross(k_s_n,Js));
    part2_2 = 1/Z0*cross(k_s_n,Jm);
    you1 = part1*(part2_1+part2_2);
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
rcs_TH = 4*pi*r^2*(abs(sum(E.*e_r_thta_s,2).^2))/Z0^2;  %36*pi*
rcs_PH = 4*pi*r^2*(abs(sum(E.*e_r_fire_s,2).^2))/Z0^2;
RCS_TH = 10*log10(rcs_TH);
RCS_PH = 10*log10(rcs_PH);
%%
%绘图
figure(1)
plot(thta_s,RCS_TH)
hold on
figure(2)
plot(thta_s,RCS_PH)


