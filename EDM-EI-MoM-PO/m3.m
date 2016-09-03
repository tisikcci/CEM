%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                          初始电压矩阵V0                                 %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear;
clc;
load('MOM.mat');
load('EH.mat');
load('Z_MOM_1.mat');
%%
%电压矩阵和外加场的设置有关
%设置外加场：
thi = 35;
phi = 0;
phi = phi*pi/180;
thi = thi*pi/180;
k_i = -[sin(thi)*cos(phi) sin(thi)*sin(phi) cos(thi)];        %入射波方向
e_i_theta = [cos(phi)*cos(thi) cos(thi)*sin(phi) -sin(thi)];  %thta极化入射
e_i_phi = [-sin(phi) cos(phi) 0];                             %fire极化入射
Pol = e_i_theta;                                              %极化方向选择
kv=k*k_i;                                                     %入射波矢
ph = 0;
phs = ph*pi/180;      %ph
%%
V0 = zeros(Edg_MOM_Total,1);                            
tic;
for m=1:Edg_MOM_Total  
    ScalarProduct=kv*Center_MOM_Plus(m,:)';             %对所有正三角面元取得编号、取得三角面元图心坐标乘以波矢并求和：即（k*r） 
    EmPlus =Pol*exp(-1j*ScalarProduct);                 %完整相位表示；得到入射场表达
    ScalarProduct=kv*Center_MOM_Minus(m,:)';   
    EmMinus=Pol*exp(-1j*ScalarProduct);     
    ScalarPlus =EmPlus* RHO_MOM_Plus(m,:)';             %E*rou     
    ScalarMinus=EmMinus*RHO_MOM_Minus(m,:)';  
    V0(m,1)=-Ed_MOM_Length(m,1)*(ScalarPlus/2+ScalarMinus/2);      
end
disp(['电压矩阵装填时间：',num2str(toc),'s']);
clear ScalarProduct EmPlus ScalarProduct EmMinus ScalarPlus ScalarMinus m
%%
tic;
%求解矩阵方程组
I0=Z_MOM_1\V0;                                                 %得到电流的表示，以公共边数目装载
disp(['矩阵求逆时间：',num2str(toc),'s']);
%%
FileName='V0_I0.mat';
save(FileName,'V0','I0');
FileName2='E_i.mat';
save(FileName2,'thi','phi','k_i','e_i_theta','e_i_phi','Pol','kv','ph','phs');