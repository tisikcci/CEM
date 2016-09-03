clear,clc
%%%%%%%%%:入射角度为（0，0）x轴极化
f=2;a=1;theta=linspace(0,2*pi,361);phi=0.5*pi;
%x=(2*pi/0.3).*f.*a;
N=length(theta);
x=(180/pi)*theta;x=x-180;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for i=1:N
    y(i)=Mie(f,a,theta(i),phi);
end
y1=pi*(a)^2;y1=10*log10(y1);y1=y1*ones(1,N);
plot(x,y,'r',x,y1,'b')
grid on