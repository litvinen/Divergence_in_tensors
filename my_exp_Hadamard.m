% Compute/approximate Hadamard_matrix-exp-function 
% Taylor: exp(x) = 1+x+x^2/2!+...x^n/n!+...
% x=A; series above is always converges, exp() is well defined function 


clc;
clear;

%A = [0.1 0.2 ; 0.3 0.4];
A = [0.5 0.2 0.3; 0.03 0.4 0.05; 0.06 0.07 0.08];

A=A*A';
E = eig(A);
rho=max(abs(E))

n=size(A,2)
I= eye(n);

B=A;%+I;
expB = exp(B)

%expmB= expm(B)

% 
% A2=A*A;
% A3=A2*A;
% A4=A3*A;
% A5=A4*A;
% A6=A5*A;
% A7=A6*A;
% TaylorB1 = I + A 
% TaylorB2 = TaylorB1 + A2/factorial(2) 
% TaylorB3 = TaylorB2 + A3/factorial(3)
% TaylorB4 = TaylorB3 + A4/factorial(4)
% TaylorB5 = TaylorB4 + A5/factorial(5) 
% TaylorB6 = TaylorB5 + A6/factorial(6) 
% TaylorB7 = TaylorB6 + A7/factorial(7)
% 
% e1=norm(expmB - TaylorB1)/norm(expmB)
% e2=norm(expmB - TaylorB2)/norm(expmB)
% e3=norm(expmB - TaylorB3)/norm(expmB)
% e4=norm(expmB - TaylorB4)/norm(expmB)
% e5=norm(expmB - TaylorB5)/norm(expmB)
% e6=norm(expmB - TaylorB6)/norm(expmB)
% e7=norm(expmB - TaylorB7)/norm(expmB)
e = ones(n,1);
Id = e*e';
A2=A.*A;
A3=A2.*A;
A4=A3.*A;
A5=A4.*A;
A6=A5.*A;
A7=A6.*A;
TaylorB1 = Id + A 
TaylorB2 = TaylorB1 + A2/factorial(2) 
TaylorB3 = TaylorB2 + A3/factorial(3)
TaylorB4 = TaylorB3 + A4/factorial(4)
TaylorB5 = TaylorB4 + A5/factorial(5) 
TaylorB6 = TaylorB5 + A6/factorial(6) 
TaylorB7 = TaylorB6 + A7/factorial(7)

e1=norm(expB - TaylorB1)/norm(expB)
e2=norm(expB - TaylorB2)/norm(expB)
e3=norm(expB - TaylorB3)/norm(expB)
e4=norm(expB - TaylorB4)/norm(expB)
e5=norm(expB - TaylorB5)/norm(expB)
e6=norm(expB - TaylorB6)/norm(expB)
e7=norm(expB - TaylorB7)/norm(expB)

