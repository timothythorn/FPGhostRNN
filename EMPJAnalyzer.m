coeff=pca(xlog');
pcavec=zeros(5000,N);
for j=1:N
for i=1:5000
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
plot(1:5000,pcavec(1:5000,1))
plot(1:5000,pcavec(1:5000,2))
plot(1:5000,pcavec(1:5000,3))

