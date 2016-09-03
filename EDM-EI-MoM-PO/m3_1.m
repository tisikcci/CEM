%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                            电流可视化程序                               %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%此程序用于在已知MOM区域电流已知的情况下进行电流可视化
%
clear;
clc;
%加载数据
load('MOM.mat');
load('V0_I0.mat');
%获得每个三角面元的电流密度
M =I0.*Ed_MOM_Length;     %电流乘公共边  Edg_MOM_Total*1
rho_A_P = RHO_MOM_Plus./repmat(2*Area_MOM(Tri_MOM_Plus,1),1,3);  %rho+/2A+  Edg_MOM_Total*3
rho_A_M = RHO_MOM_Minus./repmat(2*Area_MOM(Tri_MOM_Minus,1),1,3);%rho-/2A-
for n=1:T_NUM_MOM              %对所有面片进行遍历
    i=[0 0 0];          
    for m=1:Edg_MOM_Total      
        IE=M(m,1);         
        if(Tri_MOM_Plus(m)==n)   %找到第n个面片被当作正面片对应的位置
            i=i+IE*rho_A_P(m,:); %面片n上因被当作正面片导致的电流积累
        end
        if(Tri_MOM_Minus(m)==n)  %面片n上因被当作负面片导致的电流积累
            i=i+IE*rho_A_M(m,:);
        end
    end
    CurrentNorm(n,1)=abs(norm(i)); %面片n积累的电流累积后求模
end
toll_J = sum(CurrentNorm);%所有面片的电流累积
average_J = toll_J/T_NUM_MOM;   %求所有面片的电流平均
Jmax=max(CurrentNorm);          %找到所有面片中模值最大的电流幅值
MaxCurrent=strcat(num2str(Jmax),'[A/m]');%拼接字符串用以显示最大电流幅值
CurrentNorm1=CurrentNorm/max(CurrentNorm);%归一化
for m=1:T_NUM_MOM               %绘图
    N = t_MOM(m,:)';            %取出第m个三角面片对应的三个节点编号
    X(1:3,m)=[p_MOM(N,1)];      %三个节点的x坐标
    Y(1:3,m)=[p_MOM(N,2)];      %y
    Q(1:3,m)=[p_MOM(N,3)];      %z
end
C=repmat(CurrentNorm1,1,3)';
h=fill3(X, Y, Q, C); %绘图
colormap gray;
axis('equal');
rotate3d
FileName='current_i.mat';
save(FileName, 'CurrentNorm')