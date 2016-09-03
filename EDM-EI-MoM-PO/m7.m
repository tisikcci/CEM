%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                            MOM_PO迭代求解                                %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear;
clc;
load('Z_MOM_1.mat');
load('V0_I0.mat');
load('Matrix_PO_MOM.mat');
load('Matrix_MOM_PO.mat');
load('I_inc_k.mat');
%%
%设置迭代阀值
afa = 1e-5;  %阀值
eff = 1;     %赋迭代处置初值
I_MOM_i1 = I0;   
n=0;         %迭代次数
tic;
while(eff>afa)
    I_PO_i=Matrix_PO_MOM*I_MOM_i1+I_inc_k;     %第i次迭代PO区域的电流
    detV_i = Matrix_MOM_PO*I_PO_i+V0;          %第i次迭代MOM的修正电压
    I_MOM_i2=Z_MOM_1\detV_i;                       %第i次迭代MOM的电流系数
    det_I = I_MOM_i2-I_MOM_i1;                 %迭代前后两次电流差
    eff = norm(det_I,2)/norm(I_MOM_i1,2);      %
    I_MOM_i1 = I_MOM_i2;
    n=n+1;
    emsno(n,1)=eff;
end
disp(['迭代求解时间：',num2str(toc),'s']);
I_MOM = I_MOM_i1;
I_PO = I_PO_i;
FileName='MOMandPO_i.mat';
save(FileName,'I_MOM','I_PO');