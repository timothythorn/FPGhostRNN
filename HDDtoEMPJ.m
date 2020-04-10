syms x;
order=2;%usually 2, will use the first order-1 derivatives
epsilon=0.1%[0.1 0.05 0.01 0.2 0.5];%Epsilon beneath x axis for slowpoint node
m=order%Highest derivative used- 2 or 3 depending on if we want to specify ypp
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



