%------------------根据节点和面元信息得到几何关系---------------------%
clear all;
clc;
p0 = load('NODE.txt');
t0 = load('FACE.txt');
p0(:,1)=[];   %删除第一列
t0(:,1)=[];   %删除第一列
t = t0;
p = p0;
TrianglesTotal=length(t);          %总的三角形面元个数
clear p0 t0
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
%找到所有的公共边元，正负三角面元对
n=0;
Edge_ = [];
for m=1:TrianglesTotal             %遍历所有三角面元
   N=t(m,:)';
   for k=m+1:TrianglesTotal 
       M=t(k,:)';
       a=1-all([N-M(1) N-M(2) N-M(3)]);
       if(sum(a)==2)
           n=n+1;
           Edge_=[Edge_ M(find(a))];
           TrianglePlus(n,1)=m;
           TriangleMinus(n,1)=k;
       end
           
   end
end
EdgesTotal = length(Edge_);                  %获得公共边元数目
Edge_= Edge_';
clear n N M k a m
%公共边对应自由节点和对应法基元
A0 = t(TrianglePlus,:);
A1 = (A0 - repmat(Edge_(:,1),1,3)==0);
A2 = (A0 - repmat(Edge_(:,2),1,3)==0);
A3 = (A1+A2==0);
A4 = (A0.*A3)';
A4(A4==0)=[];
freeNodePlus = A4';
A0 = t(TriangleMinus,:);
A1 = (A0 - repmat(Edge_(:,1),1,3)==0);
A2 = (A0 - repmat(Edge_(:,2),1,3)==0);
A3 = (A1+A2==0);
A4 = (A0.*A3)';
A4(A4==0)=[];
freeNodeMinus= A4';
clear A0 A1 A2 A3 A4
%公共边长和
RHO_Plus = Center(TrianglePlus,:)-p(freeNodePlus,:);
RHO_Minus = -Center(TriangleMinus,:)+p(freeNodeMinus,:);
P1_P1 = p(Edge_(:,1),:)-p(Edge_(:,2),:);
dolp_r0 = (Center(TrianglePlus,:)+Center(TriangleMinus,:))/2;  %等效偶极子的位置
EdgeLength = sqrt(sum(P1_P1.^2,2));
clear P1_P1 p1_p2 taoA taoB taoA1 taoB1
%%
%外加场电磁参数设置
f = 3e8;                %入射波频率设置  
epsilon_ = 8.854e-012;  %自由空间介电常数
mu_ = 1.257e-006;       %自由空间磁导率
c_=1/sqrt(epsilon_*mu_);%自由空间光速
lambda =c_/f;           %入射波波长
Z0=sqrt(mu_/epsilon_);  %本征阻抗
k0 = 2*pi/lambda;       %空间自由波数
R_R = 2e3*lambda;       %空间场点距离（自行设置，但要足够远）    
%%
%利用等效偶极子方法进行散射体表面电场等效
m_m = repmat(EdgeLength,1,3).*(Center(TriangleMinus,:)-Center(TrianglePlus,:));  %电偶极矩
%for i = 1:EdgesTotal
    %r_z = repmat(Center(TrianglePlus(i,1),:),EdgesTotal,1)-dolp_r0;  %所有电偶极子辐射到正负面元上的距离矢量
    %r_f = repmat(Center(TriangleMinus(i,1),:),EdgesTotal,1)-dolp_r0;
    %r_z_2 = sum(r_z.^2,2);                                  %r_n+的平方
    %r_f_2 = sum(r_f.^2,2);                                  %r_n-的平方      
    %r_z_M = repmat(sum(r_z.*m_m,2)./r_z_2,1,3).*r_z;        %M+
    %r_f_M = repmat(sum(r_f.*m_m,2)./r_f_2,1,3).*r_f;        %M-
    %C_C_z = repmat(1./r_z_2.*(1+1./(1j*k0*sqrt(r_z_2))),1,3);
    %C_C_f = repmat(1./r_f_2.*(1+1./(1j*k0*sqrt(r_f_2))),1,3);
    %E_r_z = Z0/(4*pi)*((r_z_M-m_m).*(repmat(1j*k0./sqrt(r_z_2),1,3)+C_C_z)+2*r_z_M.*C_C_z).*repmat(exp(-1j*k0*sqrt(r_z_2)),1,3);
    %E_r_f = Z0/(4*pi)*((r_f_M-m_m).*(repmat(1j*k0./sqrt(r_f_2),1,3)+C_C_f)+2*r_f_M.*C_C_f).*repmat(exp(-1j*k0*sqrt(r_f_2)),1,3);
    %RHO_E_z = E_r_z*RHO_Plus(i,:)';
    %RHO_E_f = E_r_f*RHO_Minus(i,:)';
    %Zmn_EJ(:,i) = (EdgeLength(i,1)/2*(RHO_E_z+RHO_E_f));
    %判断各个等效偶极子之间的距离关系
    %D_d= repmat(dolp_r0(i,:),EdgesTotal,1)-dolp_r0;
    %D_D(i,:) = sqrt(sum(D_d.^2,2))';
%end
lr = 0.01*lambda;
for i = 1:EdgesTotal
    R_S = repmat(dolp_r0(i,:),EdgesTotal,1)-dolp_r0;  %两个偶极矩之间的位置关系向量表示
    r_R = sqrt(sum(R_S.^2));        %两个偶极矩之间的距离
    KL = (r_R<=lr);                 %找到两个偶极矩的距离小于阀值的部分
    r_mo = r_R+KL;                  %对其+1进行计算修正（这部分会被九点法得到的矩阵矩阵替换掉，不影响计算结果）此时这个两才表示真正的距离
    r_mn = R_S./repmat(r_mo,1,3);   %R-/R单位化矢量
    Ccc = 1./r_R.^2.*(1+1./(1j*k0*r_R));
    Zmn_EJ = Z0*exp(-1j*k0*r_mo)/(4*pi).*(m_m*m_m(i,:)'.*(1j*k0./r_mo+Ccc)-sum(m_m.*r_mn,2).*(r_mn*m_m(i,:)'.*(1j*k0./r_mo+3*Ccc)));
end
load('ZZ.mat');
ZZ = Z;
Z = Zmn_EJ';
D_D_r = (D_D>0.5*lambda);
Z = conj(Z).*D_D_r;
ZZ = ZZ.*(1-D_D_r);
Z = ZZ+Z;
%%
%入射场设置
thi0 = 35;
phi0 = 20;
phi = phi0*pi/180;
thi = thi0*pi/180;
d= -[sin(thi)*cos(phi) sin(thi)*sin(phi) cos(thi)];          %入射波方向
e_i_theta = [cos(phi)*cos(thi) cos(thi)*sin(phi) -sin(thi)];  %thta极化入射
e_i_phi = [-sin(phi) cos(phi) 0];                             %fire极化入射
Pol = e_i_theta;                                              %极化方向选择
%接收场设置
kv =k0*d;
ScalarProduct=Center(TrianglePlus,:)*kv';    
EmPlus =exp(-1j*ScalarProduct)*Pol;                
ScalarProduct=Center(TriangleMinus,:)*kv';   
EmMinus=exp(-1j*ScalarProduct)*Pol;    
ScalarPlus =sum(EmPlus.* RHO_Plus,2);            
ScalarMinus=sum(EmMinus.*RHO_Minus,2);            
V=-EdgeLength.*(ScalarPlus/2+ScalarMinus/2);
tic;
%%
%求解矩阵方程组
I=Z\V;                                                 %得到电流的表示，以公共边数目装载

%%
a = 90;
AA = 2*a+1;
thta_s = -a:a;
deg_thta_s = thta_s'*pi/180;
deg_fire_s = phi;
e_r_thta_s = [cos(deg_thta_s)*cos(deg_fire_s) cos(deg_thta_s)*sin(deg_fire_s) -sin(deg_thta_s)];%thta极化接收
e_r_fire_s = [-sin(deg_fire_s) cos(deg_fire_s) 0];  %fire极化接收
k_s = [sin(deg_thta_s)*cos(deg_fire_s) sin(deg_thta_s)*sin(deg_fire_s) cos(deg_thta_s)];%散射波方向矢量
RR = R_R*k_s;           %接收场位置
JKE = 1j*k0*Z0/(4*pi);
I_m_m = repmat(I,1,3).*m_m;
for ii=1:AA
    r = repmat(RR(ii,:),EdgesTotal,1)-dolp_r0;
    E0 = JKE./R_R.*exp(-1j*k0*(r*k_s(ii,:)'));  
    E = repmat(E0,1,3).*(repmat(sum(r.*I_m_m,2),1,3).*r./repmat(sum(abs(r).^2,2),1,3)-I_m_m);
    sum_E(ii,:) = sum(E);
end
RCS_th = 4*pi*R_R.^2*(abs(sum(sum_E.*e_r_thta_s,2)).^2);          
RCS_ph = 4*pi*R_R.^2*(abs(sum(sum_E.*repmat(e_r_fire_s,AA,1),2)).^2);
rcs_th = 10*log10(RCS_th);
rcs_ph = 10*log10(RCS_ph); 
clear E_L_i e_i_PHI sita sita_ii ii deg_thta RCS_th RCS_ph r E0 mm
%%
figure(1)
plot(-a:a,rcs_th)
hold on
figure(2)
plot(-a:a,rcs_ph)