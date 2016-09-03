%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                        装填MOM-PO区域的耦合矩阵                          %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear;
clc;
load('MOM.mat');
load('PO.mat');
load('EH.mat');
%%
%装填H_MOM_PO:即PO作用于MOM区域(散射电场)
%目的是为了获得电压矩阵det_V
Constant_part1 = -eta_/(4*pi);
Matrix_MOM_PO = zeros(Edg_MOM_Total,Edg_PO_Total);%M*K
%n_i_MOM_Plus = n_i_MOM(Tri_MOM_Plus,:);  %MOM区域正面片的法向量
%n_i_MOM_Minus = n_i_MOM(Tri_MOM_Minus,:);  %MOM区域正面片的法向量
%n_i_MOM_PM = (n_i_MOM_Plus+n_i_MOM_Minus)./repmat(sqrt(sum((n_i_MOM_Plus+n_i_MOM_Minus).^2,2)),1,3);  %面元对的法向量
tic;
for i = 1:Edg_PO_Total
    m_k =mn_PO(i,:)';   %PO区域第i个电偶极矩
    r_MOM_PO_r_n = repmat(dolp_PO_r0(i,:),Edg_MOM_Total,1); %PO区域电偶极矩的位置 M*3
    r_MOM_PO_r_n_i = dolp_MOM_r0-r_MOM_PO_r_n;   %PO区域第i个电偶极矩与MOM区域所有电偶极矩的位置关系 M*3
    norm_r_MOM_PO_r_n_i = sqrt(sum(r_MOM_PO_r_n_i.^2,2)); %PO区域第i个偶极子的位置与MOM区域所有偶极子的距离 M*1
    r_PO = r_MOM_PO_r_n_i./repmat(norm_r_MOM_PO_r_n_i,1,3);     %所有电偶极矩到第m个电偶极矩的方向矢量 R_M*3
    jk_r = exp(-1j*k*norm_r_MOM_PO_r_n_i);
    CC = 1./norm_r_MOM_PO_r_n_i.^2.*(1+1./(1j*k*norm_r_MOM_PO_r_n_i));
    Matrix_MOM_PO(:,i) = Constant_part1*jk_r.*((m_n*m_k).*(1j*k./norm_r_MOM_PO_r_n_i+CC)-...%M*1
                         sum(m_n.*r_PO,2).*(r_PO*m_k).*(1j*k./norm_r_MOM_PO_r_n_i+3*CC));
end
disp(['PO作用于MOM区域的矩阵装填时间：',num2str(toc),'s']);
FileName1 = 'Matrix_MOM_PO.mat';
save(FileName1,'Matrix_MOM_PO','-v7.3');