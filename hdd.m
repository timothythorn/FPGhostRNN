function [yi, p, pval] = hdd(x, y, yp, ypp, xi,order)
%Computes hermite interpolating polynomial coefficients via newton DD 
if order == 1
   yp = [];
   ypp = [];
   pval = length(x)*(0+1) - 1;
elseif order == 2
   ypp = [];
   pval = length(x)*(1+1) - 1;
elseif order ==3
   pval = length(x)*(2+1) - 1;
else
   fprintf('Order not valid') 
end
x = reshape(repmat(x,order,1),1,length(x)*order);
y = reshape(repmat(y,order,1),1,length(y)*order);
yp = reshape(repmat(yp,order,1),1,length(yp)*order);
ypp = reshape(repmat(ypp,order,1),1,length(ypp)*order);
n = length(x)-1;
D = zeros(n+1,n+1);
D(:,1) = y(:);
tol = n * max(x) * eps(class(x));
   
if (order == 1)
   for i = 1:n
      for j = 1:i
         D(i+1,j+1) = (D(i+1,j)-D(i,j))/(x(i+1)-x(i-j+1));
      end
   end
elseif (order == 2)
   for i = 1:n
      for j = 1:i
         h = (x(i+1)-x(i-j+1));
         if (j == 1 && h < tol)
            D(i+1,j+1) = yp(i)/1;
         else
            D(i+1,j+1) = (D(i+1,j)-D(i,j))/h;
         end
      end
   end
else
   for i = 1:n
      for j = 1:i
         h = (x(i+1)-x(i-j+1));
         if (j == 1 && h < tol)
            D(i+1,j+1) = yp(i)/1.0;
         elseif (j == 2 && h < tol)
            D(i+1,j+1) = ypp(i)/2.0;
         else
            D(i+1,j+1) = (D(i+1,j)-D(i,j))/h;
         end
      end
   end
end
p = diag(D);
for k = 1:length(xi)
   % Evaluate polynomial with the difference coefficients
   yi(k) = 0;
   for i = 1:length(p)
      term = p(i);
      for j = 1:i
         if (j == 1)
            term = term*1;
         else
            term = term*(xi(k)-x(j-1));
         end
      end
      yi(k) = yi(k) + term;
   end
end
end 