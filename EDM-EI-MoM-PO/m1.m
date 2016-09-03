%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                          网格数据处理                                   %                      
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear;
clc;
p0 = load('NODE.txt');
t0 = load('FACE.txt');
p0(:,1)=[];   %删除第一列
t0(:,1)=[];   %删除第一列
t = t0;
p = p0;
P_NUM_MOM = 831;  %MOM区域节点数
T_NUM_MOM = 1658;  %MOM区域面元数
TrianglesTotal=length(t);  %总的三角形面元个数

P_NUM_PO = 2995;   %PO区域节点数
T_NUM_PO = 5821;   %PO区域面元数
clear p0 t0
t_MOM = t(1+T_NUM_PO:end,:)-P_NUM_PO; %MOM区域面元编号
p_MOM = p(1+P_NUM_PO:end,:); %MOM区域节点编号
t_PO  = t(1:T_NUM_PO,:);  %PO区域面元编号
p_PO  = p(1:P_NUM_PO,:);   %PO区域节点编号
clear p0 t0
%%
tic;
%MOM区域每个面元的面积、图心、外法矢
A_MOM = p_MOM(t_MOM(:,1),:);                   %每个三角面元的节点坐标
B_MOM = p_MOM(t_MOM(:,2),:);
C_MOM = p_MOM(t_MOM(:,3),:);
AB_MOM= B_MOM-A_MOM;                           %三角面元边矢
BC_MOM = C_MOM-B_MOM;
CA_MOM = A_MOM-C_MOM;
ABxBC = cross(AB_MOM,BC_MOM);
Area_MOM = 0.5*sqrt(sum(ABxBC.^2,2));          %每个面元的面积
Center_MOM = 1/3*(A_MOM+B_MOM+C_MOM);                         %第m个三角面元的几何中心坐标(图心)
n_i_MOM = ABxBC./repmat(sqrt(sum(ABxBC.^2,2)),1,3);           %面元法矢
ToT_S_s_MOM = sum(Area_MOM);                                  %所有三角面元求和
clear ABxBC;
%%
%找到MOM所有的公共边元，正负三角面元对
n=0;
Edge_MOM = [];
for m=1:T_NUM_MOM                                 %遍历所有三角面元
   N=t_MOM(m,:)';
   for k=m+1:T_NUM_MOM 
       M=t_MOM(k,:)';
       a=1-all([N-M(1) N-M(2) N-M(3)]);
       if(sum(a)==2)
           n=n+1;
           Edge_MOM=[Edge_MOM M(find(a))];
           Tri_MOM_Plus(n,1)=m;
           Tri_MOM_Minus(n,1)=k;
       end
           
   end
end
Edg_MOM_Total = length(Edge_MOM);                  %获得公共边元数目
Edge_MOM= Edge_MOM';
%MOM区域公共边对应自由节点和对应法基元
A0 = t_MOM(Tri_MOM_Plus,:);
A1 = (A0 - repmat(Edge_MOM(:,1),1,3)==0);
A2 = (A0 - repmat(Edge_MOM(:,2),1,3)==0);
A3 = (A1+A2==0);
A4 = (A0.*A3)';
A4(A4==0)=[];
fr_N_MOM_Plus = A4';
A0 = t_MOM(Tri_MOM_Minus,:);
A1 = (A0 - repmat(Edge_MOM(:,1),1,3)==0);
A2 = (A0 - repmat(Edge_MOM(:,2),1,3)==0);
A3 = (A1+A2==0);
A4 = (A0.*A3)';
A4(A4==0)=[];
fr_N_MOM_Minus= A4';
%MOM区域公共边长和正负面片的图心
RHO_MOM_Plus = Center_MOM(Tri_MOM_Plus,:)-p_MOM(fr_N_MOM_Plus,:);
RHO_MOM_Minus = -Center_MOM(Tri_MOM_Minus,:)+p_MOM(fr_N_MOM_Minus,:);
PO_Ln = p_MOM(Edge_MOM(:,1),:)-p_MOM(Edge_MOM(:,2),:);   %MOM区域公共边的两个节点向量
dolp_MOM_r0 = (Center_MOM(Tri_MOM_Plus,:)+Center_MOM(Tri_MOM_Minus,:))/2;  %等效偶极子中心的位置
Ed_MOM_Length = sqrt(sum(PO_Ln.^2,2));                       %公共边长度
m_n = repmat(Ed_MOM_Length/2,1,3).*(RHO_MOM_Plus+RHO_MOM_Minus);         %MOM区域所有电偶极矩
Center_MOM_Plus = Center_MOM(Tri_MOM_Plus,:);            %正面片图心位置
Center_MOM_Minus = Center_MOM(Tri_MOM_Minus,:);          %负面片图心位置
%%
%采用九点法
%找到九个小三角形的图心：
C1=A_MOM+(1/3)*AB_MOM;      %第m个三角面元内9个子三角形的位置矢量
C2=A_MOM+(2/3)*AB_MOM;
C3=B_MOM+(1/3)*BC_MOM;
C4=B_MOM+(2/3)*BC_MOM;
C5=A_MOM-(1/3)*CA_MOM;
C6=A_MOM-(2/3)*CA_MOM;
a1=1/3*(C1+C5+A_MOM);
a2=1/3*(C1+C2+Center_MOM);
a3=1/3*(C2+C3+B_MOM);
a4=1/3*(C2+C3+Center_MOM);
a5=1/3*(C3+C4+Center_MOM);
a6=1/3*(C1+C5+Center_MOM);
a7=1/3*(C5+C6+Center_MOM);
a8=1/3*(C4+C6+Center_MOM);
a9=1/3*(C4+C6+C_MOM);
Center_=[a1 a2 a3 a4 a5 a6 a7 a8 a9];   %T_NUM_MOM*27
Center__Plus_n = Center_(Tri_MOM_Plus,:);      %所有正面片的九个图心
Center__Minus_n = Center_(Tri_MOM_Minus,:);    %所有负面片的九个图心
RHO_MOM__Plus_n= Center__Plus_n-repmat(p_MOM(fr_N_MOM_Plus,:),1,9);       %所有正面片上九个rho+
RHO_MOM__Minus_n=-Center__Minus_n+repmat(p_MOM(fr_N_MOM_Minus,:),1,9);    %所有负面片上九个rho-
disp(['MOM区域网格数据处理时间：',num2str(toc),'s']);
%%
%MOM区域数据总结：
FileName1 = 'MOM.mat';
save(FileName1,'t_MOM','p_MOM','A_MOM','B_MOM','C_MOM','AB_MOM','BC_MOM','CA_MOM',...
     'Area_MOM','Center_MOM','n_i_MOM','ToT_S_s_MOM','Tri_MOM_Plus','Tri_MOM_Minus',...
     'Edg_MOM_Total','Edge_MOM','fr_N_MOM_Plus','fr_N_MOM_Minus','dolp_MOM_r0',...
     'Ed_MOM_Length','RHO_MOM_Plus','RHO_MOM_Minus','Center_MOM_Plus','Center__Plus_n',...
     'Center_MOM_Minus','Center_','RHO_MOM__Plus_n','RHO_MOM__Minus_n','Center__Minus_n',...
     'T_NUM_MOM','P_NUM_MOM','m_n');
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                  PO区域                                 %                                                           
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tic;
%PO区域每个面元的面积、图心、外法矢
A_PO = p_PO(t_PO(:,1),:);                   %每个三角面元的节点坐标
B_PO = p_PO(t_PO(:,2),:);
C_PO = p_PO(t_PO(:,3),:);
AB_PO = B_PO-A_PO;                          %三角面元边矢
BC_PO = C_PO-B_PO;
CA_PO = A_PO-C_PO;
ABxBC = cross(AB_PO,BC_PO);
Area_PO = 0.5*sqrt(sum(ABxBC.^2,2));        %每个面元的面积
Center_PO = 1/3*(A_PO+B_PO+C_PO);                        %第m个三角面元的几何中心坐标(图心)
n_i_PO = ABxBC./repmat(sqrt(sum(ABxBC.^2,2)),1,3);       %面元法矢
ToT_S_s_PO = sum(Area_PO);                               %所有三角面元求和
r_n_PO = [AB_PO,BC_PO,CA_PO];                            %一行九列的三边矢量矩阵      
r_c_PO = [(A_PO+B_PO)/2,(B_PO+C_PO)/2,(A_PO+C_PO)/2];    %一行九列的三边中点坐标矩阵
clear ABxBC;
%%
%找到PO所有的公共边元，正负三角面元对
n=0;
Edge_PO = [];
for m=1:T_NUM_PO                                %遍历所有三角面元
   N=t_PO(m,:)';
   for k=m+1:T_NUM_PO
       M=t_PO(k,:)';
       a=1-all([N-M(1) N-M(2) N-M(3)]);
       if(sum(a)==2)
           n=n+1;
           Edge_PO=[Edge_PO M(find(a))];
           Tri_PO_Plus(n,1)=m;
           Tri_PO_Minus(n,1)=k;
       end
           
   end
end
Edg_PO_Total = length(Edge_PO);                  %获得公共边元数目
Edge_PO= Edge_PO';
clear n N M k a m
%MOM区域公共边对应自由节点和对应法基元
A0 = t_PO(Tri_PO_Plus,:);
A1 = (A0 - repmat(Edge_PO(:,1),1,3)==0);
A2 = (A0 - repmat(Edge_PO(:,2),1,3)==0);
A3 = (A1+A2==0);
A4 = (A0.*A3)';
A4(A4==0)=[]; 
fr_N_PO_Plus = A4';  
A0 = t_PO(Tri_PO_Minus,:);               %正自由节点编号
A1 = (A0 - repmat(Edge_PO(:,1),1,3)==0);
A2 = (A0 - repmat(Edge_PO(:,2),1,3)==0);
A3 = (A1+A2==0);
A4 = (A0.*A3)';
A4(A4==0)=[];
fr_N_PO_Minus= A4';                      %负自由节点编号
clear A0 A1 A2 A3 A4
%MOM区域公共边长和
Center_PO_Plus = Center_PO(Tri_PO_Plus,:);
Center_PO_Minus = Center_PO(Tri_PO_Minus,:);
RHO_PO_Plus = Center_PO_Plus-p_PO(fr_N_PO_Plus,:);          %rho+
RHO_PO_Minus = -Center_PO_Minus+p_PO(fr_N_PO_Minus,:);      %rho-
PO_Ln = p_PO(Edge_PO(:,1),:)-p_PO(Edge_PO(:,2),:);
Ed_PO_Length = sqrt(sum(PO_Ln.^2,2));                                   %ln                  
dolp_PO_r0 = (Center_PO(Tri_PO_Plus,:)+Center_PO(Tri_PO_Minus,:))/2;    %等效偶极子的位置
mn_PO = repmat(Ed_PO_Length/2,1,3).*(RHO_PO_Plus+RHO_PO_Minus);         %PO区域所有电偶极矩
disp(['PO区域网格数据处理时间：',num2str(toc),'s']);
clear P1_P1 p1_p2 taoA taoB taoA1 taoB1 PO_Ln a 
%%
%PO区域数据总结：
FileName2 = 'PO.mat';
save(FileName2,'t_PO','p_PO','A_PO','B_PO','C_PO','AB_PO','BC_PO','CA_PO',...
     'Area_PO','Center_PO','n_i_PO','ToT_S_s_PO','Tri_PO_Plus','Tri_PO_Minus',...
     'Edg_PO_Total','Edge_PO','fr_N_PO_Plus','fr_N_PO_Minus','dolp_PO_r0',...
     'Ed_PO_Length','Center_PO_Plus','Center_PO_Minus','RHO_PO_Plus','RHO_PO_Minus',...
     'P_NUM_PO','T_NUM_PO','mn_PO');