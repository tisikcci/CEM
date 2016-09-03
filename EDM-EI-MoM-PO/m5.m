%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                        装填耦合矩阵Z_PO_MOM                              %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear;
clc;
load('MOM.mat');
load('PO.mat');
load('PO_I_s.mat');
load('EH.mat');
%%
%装填H_PO_MOM:即MOM作用于PO区域(散射磁场)
%目的是为了获得MOM对PO耦合贡献的电流I_1_PO
Constant1 = 1j*k/(4*pi);
Matrix_PO_MOM = zeros(Edg_PO_Total,Edg_MOM_Total); %用来存储计算MOM作用于PO时的电流系数的矩阵 K*M
tic;
for i = 1:Edg_MOM_Total
    m_n_i =repmat(m_n(i,:),Edg_PO_Total,1);   %K*3  
    sita_k_n = sita_PO_MOM(i,:)';  %MOM区域第i个公共边对PO区域所有面片的遮蔽系数判断 P*1
    sita_k_n_P = sita_k_n(Tri_PO_Plus,1);  %PO正面片遮蔽系数    K*1
    sita_k_n_M = sita_k_n(Tri_PO_Minus,1); %PO负面片遮蔽系数    K*1
    %偶极矩
    r_PO_MOM_r_n=repmat(dolp_MOM_r0(i,:),Edg_PO_Total,1); %取出MOM区域第i个偶极子的位置  K*3  
    r_PO_MOM_r_n_p = dolp_PO_r0-r_PO_MOM_r_n;             %MOM区域第i个偶极子的位置与PO区域所有偶极子的位置关系  K*3
    norm_r_PO_MOM_r_n_p = sqrt(sum(r_PO_MOM_r_n_p.^2,2)); %MOM区域第i个偶极子的位置与PO区域所有偶极子的距离   K*1
    jk_r = exp(-1j*k*norm_r_PO_MOM_r_n_p);                %相位关系 K*1
    CC = 1./norm_r_PO_MOM_r_n_p.^2.*(1+1./(1j*k*norm_r_PO_MOM_r_n_p));  %K*1
    M_part1 = 0.5*Constant1*(sita_k_n_P+sita_k_n_M).*CC.*jk_r;
    M_part2 = sum(tao_zf.*cross(n_i_PO_PM,cross(m_n_i,r_PO_MOM_r_n_p)),2);
    Matrix_PO_MOM(:,i) = M_part1.*M_part2;
end
disp(['MOM作用于PO区域的矩阵装填时间：',num2str(toc),'s']);
FileName1 = 'Matrix_PO_MOM.mat';
save(FileName1,'Matrix_PO_MOM','-v7.3');