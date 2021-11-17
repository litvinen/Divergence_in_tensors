% Compute/approximate Hadamard matrix-inverse-function
% Taylor: 1/(1-x) = 1+x+x^2+...x^n+...
% Taylor: 1/x = 1/(1-(1-x)) = 1+(1-x)+(1-x)^2+...(1-x)^n+...
% x=A; converges if  || A-I || <1


clc;
clear;
format long;
%A = [0.1 0.2 ; 0.3 0.4];
%A = [1.1 1.002 1.003; 1.002 1.2 1.005; 1.003 1.005 1.1];

%A = A+0.1*randn(3);
%A=A*A';

d=5;
n=3;
r=4;
Id = tt_ones(n, d);
xmin=2*ones(1,n);

%A = tt_random(n,d,r); 
A = tt_x(n, d, xmin);
A = A/(n^d);
%A=A.*A/n^d;
%A=A.*B;
trunc_eps = 0.0001;
%A = round(A, eps);

%scal_prod = dot(A, Id)

% sum=0.0;
% for i1=1:n
%   for i2=1:n
%     for i3=1:n
%       for i4=1:n
%         for i5=1:n
%            sum = sum  +  A(i1,i2,i3,i4,i5);
%                
%         end
%       end
%     end
%   end
% end  
% sum
% sz=size(A)




%Check necessary condition
alpha = 1/power(norm(A),2);
alpha = 0.99*alpha;
V0 = alpha*A;
ER=Id - A.*V0;
% maxel=ER(1,1,1,1,1);
% for i1=1:n
%   for i2=1:n
%     for i3=1:n
%       for i4=1:n
%         for i5=1:n
%             %ER(i1,i2,i3,i4,i5);
%            if maxel< ER(i1,i2,i3,i4,i5)
%                maxel= ER(i1,i2,i3,i4,i5);
%            end   
%                
%         end
%       end
%     end
%   end
% end  
% maxel
err_check = norm(ER)
err_check2 = norm(Id - alpha*A.*A)



%Compute exact = 1./A via iterations
V = iteration_with_truncation(@function_psi, V0, A, Id, trunc_eps)
abs_err = norm(Id - A.*V)
abs_err_sc = norm(Id - A.*V)/n^d




%================ Let us check how cross() works
%xmin=2*ones(1,d);
%W=tt_x(n, d, xmin)/n^d;

%h=1.0./(n - 1);
%x= tt_x (2 ,d)*h;
fun_inv = @(x) 1./x;
invW_cross = funcrs2(A, fun_inv, 1e-12, A ,18);
display(invW_cross)
dif_a = norm(invW_cross - V)
dif_r = norm(invW_cross - V)/norm(V)
abs_err_cross = norm(Id - A.*invW_cross)
abs_err_iter = norm(Id - A.*V)



%+++++++++++++++++




%e2= Id - A.*V

function y = iteration_with_truncation(function_iteration, xV, xA, I, eps)
  k=0;
  while (k<100)
      xV = function_iteration(xV, xA, I);
      %add here tunctation procedure
      %err = norm(exact - xV)/norm(exact)
      xV=round(xV, eps);
      k=k+1
  end
  y = xV;

end


function Y = function_psi(V, A, I)
   Y = V.*(2*I - A.*V);
end
