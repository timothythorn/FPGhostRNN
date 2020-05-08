eps=[.02];
Is=[-9.0 -9.25 -9.5 -9.75];
plots=[];
names=[];
N=50;
x0=0*randn(N,1);%0.001*randn(N,1)
nt=20000;
%I=[0*ones(.05*nt,1);-.001*ones(.001*nt,1);0*ones(.099*nt,1);-1*ones(.1*nt,1);0*ones(.1*nt,1);1*ones(.1*nt,1);0*ones(.1*nt,1);-10*ones(.25*nt,1);0*ones(.2*nt,1)];
I=[0*ones(.05*nt,1);.001*ones(.001*nt,1);0*ones(.099*nt,1);-9.25*ones(.25*nt,1);0*ones(.6*nt,1)];%0*ones(.24*nt,1);0*ones(.35*nt,1);0*ones(.35*nt,1)];
dt=.01;
SimTime=dt*nt

for qq=1:size(Is,2)
I=[0*ones(.05*nt,1);.001*ones(.001*nt,1);0*ones(.099*nt,1);Is(qq)*ones(.25*nt,1);0*ones(.6*nt,1)];%0*ones(.24*nt,1);0*ones(.35*nt,1);0*ones(.35*nt,1)];    
syms x;
order=2;%usually 2, will use the first (order-1) derivatives
epsilon=eps(1);%Epsilon beneath x axis for slowpoint node
m=order;%Highest derivative used- 2 or 3 depending on if we want to specify ypp
xi=linspace(-1,5.5,200); %For plotting the result
xs=[0 1 2 3 4];  %Node x locations
%node y locations
yp=[-5 0 0 0 -5]; %first derivative at each node (define stability)
ypp=[1 -1 1 -10 -1];  %Specify second derivative, if using order 3

% for e=1:size(epsilon,2)
ys=[0 1 -1 0-epsilon -2.5];
nn=size(xs,2);
for j=1:nn     %Lookup table for crafting the symbolic polynomial from the result
    for k=1:m
        xnm(k+(j-1)*m)=xs(j);
    end
end
[yi, P, Pv] = hdd(xs,ys,yp,ypp,xi,order); %Call to Hermite Divided Difference function to obtain polynomial

f=P(1);  %Translate CoEffs into symbolic expression f
for i=2:length(P)
    for j=1:i
        if j==1
            term=1;
        else
            term=term.*(x-xnm(j-1));
        end
    end
f=f+P(i)*term;
end

fexact=expand(f)
fexp=vpa(fexact)
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
A=reshape(A,N^2,1);
B=reshape(B,N^2,1);
A=diag(A);B=diag(B);
W=linsolve(B,A);
W=diag(W);
W=reshape(W,N,N);



clear xlog
xlog(:,1)=x0;%zeros(N,1);%
kmod=1 ;%Plot every kth cycle
stepscale=1 ;%multiply W by this each step
W=W;
%I=[0*ones(.15*nt,1);-10*ones(.25*nt,1);0*ones(.2*nt,1);-20*ones(.2*nt,1);0*ones(.2*nt,1)];


for i=2:nt
    xlog(:,i)=xlog(:,i-1)+dt*(-xlog(:,i-1)+W'*tanh(xlog(:,i-1))+I(i));
end


y=1:nt;

W=stepscale*W;



coeff=pca(xlog');
pcavec=zeros(nt,N);
for j=1:N
for i=1:nt
    pcavec(i,j)=coeff(j,:)*xlog(:,i);
end
end


figure(5)
hold on
plots(end+1)=plot(1:nt,pcavec(1:nt,1));
names{end+1}=['I = '  num2str(Is(qq))];
%plot(1:nt,pcavec(1:nt,2))
%plot(1:nt,pcavec(1:nt,3))
WLog(:,:,qq)=W;
XLogT(:,:,qq)=xlog;

end

legend(plots,names)
title('')
xlabel('t')
ylabel('PC1')