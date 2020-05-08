figure(1)
nt=20000;
dt=.005;
clear xlog
xlog(:,1)=0.001*randn(N,1);%zeros(N,1);%
kmod=1 ;%Plot every kth cycle
stepscale=1 ;%multiply W by this each step
W=W;
I=[0*ones(.15*nt,1);-10*ones(.25*nt,1);0*ones(.2*nt,1);10*ones(.2*nt,1);0*ones(.2*nt,1)];
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
pcavec=zeros(nt,N);
for j=1:N
for i=1:nt
    pcavec(i,j)=coeff(j,:)*xlog(:,i);
end
end
t1=750;
t2=2000;
t3=3000;
t4=4000;
tf=5000;
figure(4)
%plot3(pcavec(1:t1,1),pcavec(1:t1,2),pcavec(1:t1,3),'ro')

plot3(pcavec(t1:t2,1),pcavec(t1:t2,2),pcavec(t1:t2,3),'bo')
hold on
plot3(pcavec(t2:t3,1),pcavec(t2:t3,2),pcavec(t2:t3,3),'go')
plot3(pcavec(t3:t4,1),pcavec(t3:t4,2),pcavec(t3:t4,3),'ko')
plot3(pcavec(t4:tf,1),pcavec(t4:tf,2),pcavec(t4:tf,3),'mo')
xlabel('PC 1')
ylabel('PC 2')
zlabel('PC 3')

figure(5)
hold on
plot(1:nt,pcavec(1:nt,1))
plot(1:nt,pcavec(1:nt,2))
plot(1:nt,pcavec(1:nt,3))

