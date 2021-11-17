% Compute/approximate Hadamard matrix-inverse-function
% Taylor: 1/(1-x) = 1+x+x^2+...x^n+...
% Taylor: 1/x = 1/(1-(1-x)) = 1+(1-x)+(1-x)^2+...(1-x)^n+...
% x=A; converges if  || A-I || <1


clc;
clear;
format long;
%A = [0.1 0.2 ; 0.3 0.4];
A = [1.1 1.002 1.003; 1.002 1.2 1.005; 1.003 1.005 1.1];

%A = A+0.1*randn(3);
%A=A*A';

n=size(A,2)

% detA=det(A)
% detAmI=det(A-I)
% 
% E = eig(A-I);
% rho=max(abs(E))
% 2/rho
% alpha=1/rho;


e = ones(n,1);
Id = e*e';

alpha = 1/power(norm(A,Inf),2);
alpha = 0.99*alpha;
V0 = alpha*A;
err_check = norm(Id - A.*V0, Inf)
err_check2 = norm(Id - alpha*A.*A, Inf)

exact=1./A
V = iteration_with_truncation(@function_psi, V0, A, Id, exact )
e1= V-exact
e2= Id - A.*V
function y = iteration_with_truncation(function_iteration, xV, xA, I, exact)

  err=10000;
  k=0;
  while ((err > 0.001)&&(k<500))
      xV = function_iteration(xV, xA, I);
      %add here tunctation procedure
      err = norm(exact - xV)/norm(exact)
      k=k+1
  end
  y = xV;

end


function Y = function_psi(V, A, I)
   Y = V.*(2*I - A.*V);
end
