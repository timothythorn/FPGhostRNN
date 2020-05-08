N=40;
g_gg=1;
Tau=.1;
k=64;

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

RandMat=randn(N,N)*0;
A=A+RandMat;
A=reshape(A,1600,1);
B=reshape(B,1600,1);
A=diag(A);B=diag(B);
W=linsolve(B,A);
W=diag(W);
W=reshape(W,N,N);
%WBig=linsolve(B,A+RandMat);
%size(WBig)
%norm(WBig)
%Wnormed=WBig./norm(WBig)*10;
%norm(Wnormed)

%W=sqrt(g_gg^2/N)*randn(N);
%x=randn(N,1); 
%I=0;
%[U,E]=eig(W)

