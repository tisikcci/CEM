%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                           装填自阻抗矩阵                                  %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear;
clc;
%导入数据
load('MOM.mat');
%%
%电磁基本参数设置
f = 6e8;                %外加电磁波频率
epsilon_ = 8.854e-012;  
mu_ = 1.257e-006;
c_=1/sqrt(epsilon_*mu_);
eta_=sqrt(mu_/epsilon_);   %自由空间介电常数
%同期变量（用于计算阻抗矩阵）
omega       =2*pi*f;                                           
k           =omega/c_; 
K           =1j*k;
lambda = c_/f;  %波长
%%
deta = 0.15*c_/f;            %用以判断是否使用等效偶极矩方法，选取距离超过0.15个波长的
Constant_part1 = -1j*omega*mu_/(144*pi);
Constant_part2 = -4/(k^2);
%装填Z_MOM_1
Z_MOM_1 = zeros(Edg_MOM_Total,Edg_MOM_Total)+1j*zeros(Edg_MOM_Total,Edg_MOM_Total);
tic;
for m=1:Edg_MOM_Total
    rho_m_P = RHO_MOM_Plus(m,:);           %rho_m+      第m个正面片的rho_m+ 3*1
    rho_m_M = RHO_MOM_Minus(m,:);          %rho_m-      第m个负面片的rho_m-
    lm = Ed_MOM_Length(m,1);               %lm          第m个偶极矩对应的公共边的长度
    %九点法准备
    r_m_Plus_n =repmat(Center_MOM_Plus(m,:),Edg_MOM_Total,9);     %r_m+  Edg_MOM_Total*27
    r_m_Minus_n = repmat(Center_MOM_Minus(m,:),Edg_MOM_Total,9);  %r_m-  Edg_MOM_Total*27
    r_mn_PP= r_m_Plus_n-Center__Plus_n;    %r_m+-r_k+   所有其他正面片的九个图心到第m个公共边对应的正面片的图心的矢量 Edg_MOM_Total*27
    r_mn_PM = r_m_Plus_n-Center__Minus_n;  %r_m+-r_k-   所有其他负面片的九个图心到第m个公共边对应的正面片的图心的矢量 Edg_MOM_Total*27
    r_mn_MP = r_m_Minus_n-Center__Plus_n;  %r_m--r_k+   所有其他正面片的九个图心到第m个公共边对应的负面片的图心的矢量 Edg_MOM_Total*27
    r_mn_MM = r_m_Minus_n-Center__Minus_n; %r_m--r_k-   所有其他正面片的九个图心到第m个公共边对应的正面片的图心的矢量 Edg_MOM_Total*27
   %针对第m个公共边对所有个公共边的作用进行遍历：
   for n=1:Edg_MOM_Total  
       ln = Ed_MOM_Length(n,1);                %ln  
       r_mk_PP = reshape(r_mn_PP(n,:),3,9)';   %9*3 两对三角形的关系 取出第n个公共边对应的正面片九点和第m个公共边对应的正面片的向量关系
       r_mk_PM = reshape(r_mn_PM(n,:),3,9)';   %9*3
       r_mk_MP = reshape(r_mn_MP(n,:),3,9)';   %9*3
       r_mk_MM = reshape(r_mn_MM(n,:),3,9)';   %9*3
       norm_r_mk_PP = sqrt(sum(r_mk_PP.^2,2)); %9*1  对上面得到的向量取模 得到距离关系  
       norm_r_mk_PM = sqrt(sum(r_mk_PM.^2,2)); %9*1  
       norm_r_mk_MP = sqrt(sum(r_mk_MP.^2,2)); %9*1
       norm_r_mk_MM = sqrt(sum(r_mk_MM.^2,2)); %9*1
       g_PP = exp(-1j*k*norm_r_mk_PP)./norm_r_mk_PP;  %9*1 得到相位关系
       g_PM = exp(-1j*k*norm_r_mk_PM)./norm_r_mk_PM;  %9*1
       g_MP = exp(-1j*k*norm_r_mk_MP)./norm_r_mk_MP;  %9*1
       g_MM = exp(-1j*k*norm_r_mk_MM)./norm_r_mk_MM;  %9*1
       rho_k_P = reshape(RHO_MOM__Plus_n(n,:),3,9)';  %rho_k+ = 9*3  第n个公共边对应的正面片的rho_k+ 
       rho_k_M = reshape(RHO_MOM__Minus_n(n,:),3,9)'; %rho_k- = 9*3  第n个公共边对应的负面片的rho_k-
       rho_g_PP_rhom_z = sum((rho_k_P*rho_m_P').*g_PP);     %1*3
       rho_g_PM_rhom_z = sum((rho_k_M*rho_m_P').*g_PM);     %1*3
       rho_g_MP_rhom_f = sum((rho_k_P*rho_m_M').*g_MP);     %1*3
       rho_g_MM_rhom_f = sum((rho_k_M*rho_m_M').*g_MM);     %1*3

       Z_MOM_1(m,n)=Constant_part1*lm*ln*((rho_g_PP_rhom_z+rho_g_PM_rhom_z)+(rho_g_MP_rhom_f+rho_g_MM_rhom_f)+...
                    ((sum(g_PP)-sum(g_PM))-(sum(g_MP)-sum(g_MM)))*Constant_part2); %使用传统的九点法进行装填
   end
end
disp(['阻抗矩阵Z_MOM装填时间是：',num2str(toc)]);
%%
FileName1='Z_MOM_1.mat';
save(FileName1,'Z_MOM_1');
FileName2='EH.mat';
save(FileName2,'f','omega','mu_','epsilon_','c_', 'k','eta_','lambda');