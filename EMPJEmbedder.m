N=40;
g_gg=1;
Tau=.1;
%x = 0:0.5:4.5;
%y = subs(f);
k=64;
% th=0:pi/k:pi;
% ph=0:2*pi/k:2*pi;
% x=sin(th).*cos(ph);
% y=cos(th).*sin(ph);
% z=cos(th);
% X=[x;y;z];
[q,r]=gramschmidt(randn(N,1));
A=zeros(N,N);
B=zeros(N,N);
xi=0.01:(4.5)/N:4.5;
for i=1:N
    x=xi(i);
    Ut=x*q;
    y=subs(f);
    A(i,:)=(Ut*(1-tanh(x)^2))';
    B(i,:)=Tau*(y+(1/Tau)*1)*Ut;
end

RandMat=0;%randn(N,N)*0.01;
WBig=(B-RandMat)\(A+RandMat);
size(WBig)
norm(WBig)
Wnormed=WBig./norm(WBig);
norm(Wnormed)

%W=sqrt(g_gg^2/N)*randn(N);
%x=randn(N,1);
%I=0;
%[U,E]=eig(W)

