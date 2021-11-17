% Compute/approximate matrix-log-function 
%% Taylor: log(1-x) = -x-x^2/2- x^3/3 - x^4/4 +, where x = 1-A and (1-x) should be >0
%% Taylor: log(1+x) = x-x^2/2+ x^3/3 - x^4/4 +
%% Taylor: log(x)=log(1+(x-1)) = (x-1)-(x-1)^2/2+ (x-1)^3/3 - (x-1)^4/4 +
% Taylor: log(A)=log(1-x) = -(x+x^2/2+ x^3/3 + x^4/4 +...)
% x=1-A; series above converges when \rho(x)<1  
% where 1 is a matrix consisting of all "1"

% Gregory's series
clc;
clear;



%A = [0.1 0.2 ; 0.3 0.4];
W =  [0.5 0.998 0.997; 0.998 0.4 0.995; 0.993 0.995 0.3];
n=size(W,2)
e = ones(n,1);
Id = e*e';

exact = log(W)

A = Id - W
B = Id + W
invHB = 1./B
C0 = A.*invHB;
% temp = C0.*C0;
% C1 = temp.*C0;
% C2 = C1.*temp;
% C3 = C2.*temp;
% C4 = C3.*temp; %  =C0.*temp^4=C0^9

K=4
Y= iteration_with_truncation(K, C0, exact, Id)
function y = iteration_with_truncation(K, C0, exact, Id)
  sum = C0;
  for k = 1:K
    prod = Id;
    for j=1:(2*k+1)
        prod = prod .* C0;
        
    end    
   sum = sum + prod/(2*k+1);
   err = norm(exact - (-2)*sum)/norm(exact)
  end 
  y = sum;
end
% Gregory0 = C0
% Gregory1 = Gregory0 + C1/3 
% Gregory2 = Gregory1 + C2/5 
% Gregory3 = Gregory2 + C3/7
% Gregory4 = Gregory3 + C4/9 
% 
% 
% e1=norm(exact - (-2)*Gregory0)/norm(exact)
% e2=norm(exact - (-2)*Gregory1)/norm(exact)
% e3=norm(exact - (-2)*Gregory2)/norm(exact)
% e4=norm(exact - (-2)*Gregory3)/norm(exact)
% e5=norm(exact - (-2)*Gregory4)/norm(exact)


