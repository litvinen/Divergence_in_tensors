% Compute/approximate Hadamard (pointwise) log-function 
%% Taylor: log(1+x) = x-x^2/2+ x^3/3 - x^4/4 +
%% Taylor: log(x)=log(1+(x-1)) = (x-1)-(x-1)^2/2+ (x-1)^3/3 - (x-1)^4/4 +
% Taylor: log(A)=log(1-x) = -(x+x^2/2+ x^3/3 + x^4/4 +...)
% x=1-A; series above converges when \rho(x)<1  
% where 1 is a matrix consisting of all "1"

clc;
clear;



%A = [0.1 0.2 ; 0.3 0.4];
A =  [0.5 0.998 0.997; 0.998 0.4 0.995; 0.993 0.995 0.3];

%A=A*A';
E = eig(A);
rho=max(abs(E))

n=size(A,2)
e = ones(n,1);
Id = e*e';



X=A-Id
norm(X)
E=eig(X);
rho = max(abs(E))
%exact_sol 
logB= log(A)  % This is pointwise logarithm
%    0.1310   -3.6497   -3.1236
%   -3.6497    0.0050   -4.7560
%   -3.1236   -4.7560    0.0148
%logmB= logm(B)

  %  0.1299    0.0241    0.0408
  %  0.0241    0.0046    0.0080
  %  0.0408    0.0080    0.0139
    
X2=X.*X;
X3=X2.*X;
X4=X3.*X;
X5=X4.*X;
X6=X5.*X;
X7=X6.*X;
TaylorB1 = X
TaylorB2 = TaylorB1 - X2/2 
TaylorB3 = TaylorB2 + X3/3 
TaylorB4 = TaylorB3 - X4/4
TaylorB5 = TaylorB4 + X5/5 
TaylorB6 = TaylorB5 - X6/6 
TaylorB7 = TaylorB6 + X7/7 

e1=norm(logB - TaylorB1)/norm(logB)
e2=norm(logB - TaylorB2)/norm(logB)
e3=norm(logB - TaylorB3)/norm(logB)
e4=norm(logB - TaylorB4)/norm(logB)
e5=norm(logB - TaylorB5)/norm(logB)
e6=norm(logB - TaylorB6)/norm(logB)
e7=norm(logB - TaylorB7)/norm(logB)


