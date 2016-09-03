function [delta]=Mie(f,a,theta,phi)
%%%%%:求球beseel及其导数
lambda=0.3/f;K=2*pi/lambda;ka=K*a;
kr=1e4;
N=200;
sk1=j*exp(-j*ka)*cos(phi)/kr;%%%%%%:令kr=10000
sk2=-j*exp(-j*ka)*sin(phi)/kr;
n=linspace(1,N,N);
T1=(2*n+1)./(n+1)./n;T=zeros(1,N);
for i=1:N
    T(i)=(-1)^i;
end
T=T1.*T;
J=besselj(n+0.5,ka);H=besselh(n+0.5,2,ka);s=sqrt(pi/2./ka);
J=s*J;H=s*H;%%%%1*1001
J1=besselj(n-0.5,ka);H1=besselh(n-0.5,2,ka);
J1=s*J1;H1=s*H1;
an=J./H;
bn=(ka*J1-n.*J)./(ka*H1-n.*H);
%%%%%%%:求Pn(cos(theta))及其导数
x=linspace(theta-0.1,theta+0.1,3);x1=cos(x);
P=zeros(1,N);d_P=zeros(1,N);sh=eps+sin(theta);
for i=1:N
    a=legendre(i,x1);a=a(2,:);
    a1=a(2)-a(1);a2=a(3)-a(2);
    a1=a1/0.1;a2=a2/0.1;
    d_P(i)=(a1+a2)/2;
    P(i)=a(2)/sh;
end
%%%%%%%%:求E_theta、E_phi
E_theta=sk1*T.*(bn.*d_P-an.*P);
E_phi=sk2*T.*(bn.*P-an.*d_P);
E_theta=abs(sum(E_theta));
E_phi=abs(sum(E_phi));
E_theta=E_theta^2;E_phi=E_phi^2;
E=E_theta+E_phi;
delta=4*pi*(kr/K)^2*E;
delta=10*log10(delta);