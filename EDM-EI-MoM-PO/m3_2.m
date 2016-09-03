%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%          根据MOM自阻抗矩阵得到的电流系数计算MOM区域单独存在时的RCS         %        
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%加载数据：
load('EH.mat');
load('E_i.mat');
load('MOM.mat');
load('V0_I0.mat');
%初始数据赋值：
rr = 200*lambda;  %远场距离
th = 0:360;
AA = length(th);
E_s = zeros(AA,3);  %用来装填各个角度的散射接收场
EE =zeros(AA,1);
Etotal=zeros(AA,1);
ths = th'*pi/180;                                     %AA个散射接收角     AA*1
ks = [sin(ths)*cos(phs) sin(ths)*sin(phs) cos(ths)];  %散射场单位矢量     AA*3
RR = rr*ks;                                           %远场位置矢量       AA*3
er0 = [cos(ths)*cos(phs) cos(ths)*sin(phs) -sin(ths)];%theta极化方向
er1 = repmat([-sin(phs) cos(phs) 0],AA,1);            %phi极化           AA*3
er =er0;                                              %thta极化接收
mm = repmat(I0,1,3).*m_n;
Constant1 = eta_/(4*pi);   
for i=1:AA                  %遍历散射接收角
    RR_n = repmat(RR(i,:),Edg_MOM_Total,1)-dolp_MOM_r0;  %当散射接收角确定后，每个偶极矩的散射源点到场点的矢量位置
    norm_RR_n = sqrt(sum(RR_n.^2,2));                    %距离
    CC = 1./norm_RR_n.^2.*(1+1./(1j*k*norm_RR_n));
    MM = repmat(sum(RR_n.*mm,2)./norm_RR_n.^2,1,3).*RR_n;
    E = Constant1*((MM-mm).*repmat(1j*k./norm_RR_n+CC,1,3)+2*MM.*repmat(CC,1,3)).*repmat(exp(-1j*k*norm_RR_n),1,3);
    E_r = E.*repmat(norm_RR_n,1,3);
    E_s(i,:) = sum(E);                                    %同一散射接收角下的所有面片的散射场叠加
end
RCS = 4*pi*rr^2*abs(sum(E_s.*er,2).^2);
dB_RCS = 10*log10(RCS);
plot(th,dB_RCS)