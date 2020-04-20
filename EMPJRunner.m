nt=5000;
dt=.01;
clear xlog
xlog(:,1)=0.1*randn(N,1);%zeros(N,1);%
kmod=1 ;%Plot every kth cycle
stepscale=2 ;%multiply W by this each step
W=Wnormed;
I=[0*ones(.15*nt,1);0.5*ones(.25*nt,1);0*ones(.2*nt,1);-0.5*ones(.2*nt,1);0*ones(.2*nt,1)];
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