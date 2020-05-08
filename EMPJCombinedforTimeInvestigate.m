syms x;
order=2;%usually 2, will use the first (order-1) derivatives
epsilon=.1;%[0.1 0.05 0.01 0.2 0.5];%Epsilon beneath x axis for slowpoint node
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


%Simplify output, copy to clipboard
% fsimp=simplify(f,'Steps',50,'All',true);  Defunct, replaced below
% if size(fsimp,1)>1
%     clipboard('copy',strcat('simplify (',' ',char(fsimp(2)),')'))
% else
%     clipboard('copy',strcat('simplify (',' ',char(fsimp),')'))
% end
fexact=expand(f)
fexp=vpa(fexact)

clipboard('copy',strcat(char(fexact)))
fprintf('f copied to clipboard\n')
figure(1)
hold on
plot(xi,yi)
ylim([-2,4])
plot(xi,yi+2.5)
yline(0);
for i=1:size(ys,2)
    xline(xs(i));
end
hold off

I_on=2
nsteps=500;
dt=0.1;
xt(1)=0;
I(1:ceil(nsteps/4))=0;
I(ceil(nsteps/4)+1:ceil(nsteps/2))=I_on;
I(ceil(nsteps/2)+1:nsteps)=0;
for step=2:nsteps
    x=xt(step-1);
    xt(step)=xt(step-1)+dt*(subs(f)+I(step));
end


figure(2)
hold on
plot(1:nsteps,xt,'o')
hold off
%end

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

RandMat=randn(N,N)*epsilon^2;
WBig=(B-RandMat)\(A+RandMat);
size(WBig)
norm(WBig)
Wnormed=WBig./norm(WBig);
norm(Wnormed)

figure(1)
nt=5000;
dt=.01;
clear xlog
xlog(:,1)=0.1*randn(N,1);%zeros(N,1);%
kmod=1 ;%Plot every kth cycle
stepscale=2 ;%multiply W by this each step
W=WBig;%Wnormed;
I=[0*ones(.15*nt,1);0.5*ones(.25*nt,1);0*ones(.2*nt,1);1*ones(.2*nt,1);0*ones(.2*nt,1)];
for k=1:1
if mod(k,kmod)==0
    figure(k)
end
hold on

for i=2:nt
    xlog(:,i)=xlog(:,i-1)+dt*(-xlog(:,i-1)+W'*tanh(xlog(:,i-1))+I(i));
end


y=1:nt;
for j=1:N
for ii=1:10
    if mod(k,kmod)==0
    plot(y,xlog(j,:))
    nw=norm(W);
    title(['Activity for norm W=' ,string(nw)])
    end
end
end
W=stepscale*W;
end
coeff=pca(xlog');
pcavec=zeros(5000,N);
for i=1:5000
    pcavec(i,:)=coeff'*xlog(:,i);
end
t1=750;
t2=2000;
t3=3000;
t4=4000;
tf=5000;
figure(4)
plot3(pcavec(1:t1,1),pcavec(1:t1,2),pcavec(1:t1,3),'ro')
hold on
plot3(pcavec(t1:t2,1),pcavec(t1:t2,2),pcavec(t1:t2,3),'bo')
plot3(pcavec(t2:t3,1),pcavec(t2:t3,2),pcavec(t2:t3,3),'go')
plot3(pcavec(t3:t4,1),pcavec(t3:t4,2),pcavec(t3:t4,3),'ko')
plot3(pcavec(t4:tf,1),pcavec(t4:tf,2),pcavec(t4:tf,3),'mo')
xlabel('PC 1')
ylabel('PC 2')
zlabel('PC 3')

figure(5)
hold on
%plot(750:3000,pcavec(750:3000,1))
plot(750:1000,pcavec(750:1000,2))
plot(750:1000,pcavec(750:1000,3))