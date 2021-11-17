% Compute/approximate Tensor-log-function 
%% Taylor: log(1-x) = -x-x^2/2- x^3/3 - x^4/4 +, where x = 1-A and (1-x) should be >0
%% Taylor: log(1+x) = x-x^2/2+ x^3/3 - x^4/4 +
%% Taylor: log(x)=log(1+(x-1)) = (x-1)-(x-1)^2/2+ (x-1)^3/3 - (x-1)^4/4 +
% Taylor: log(A)=log(1-x) = -(x+x^2/2+ x^3/3 + x^4/4 +...)
% x=1-A; series above converges when \rho(x)<1  
% where 1 is a matrix consisting of all "1"

% Use Gregory's series to compute tensor LOG() function
clc;
clear;

d=3;
n=3;
r=2;
eps_exp=1e-3;
trunc_r= 4;
trunc_eps =0.01;
Id = tt_ones(n, d);
%expId = tt_exp(Id, eps_exp);

xmin=2*ones(1,d);
W=tt_x(n, d, xmin); % Generate a positive tensor of size n^d
rank(W) % print some ranks
display(W) % print info about tensor

%Compute exp(), this function exists
%expW = tt_exp(W, eps_exp); % Need it for testing purposes exp(log(x))= = log(exp(x)) = x
%rank(expW)

%===== Let us check how tt_exp() works
%h=1.0./(n - 1);
%x= tt_x (2 ,d)*h;
fun_exp = @(x) exp(x);
expW_cross = funcrs2(W, fun_exp, 1e-12, W ,18);  % use the cross method to compute f(w)
display(expW_cross)
fun_log = @(x) log(x);
logW_cross = funcrs2(W, fun_log, 1e-12, W ,18);  % use the cross method to compute f(w)
display(logW_cross)
W_appr = funcrs2(expW_cross, fun_log, 1e-12, expW_cross ,18);
norm(W - W_appr)
norm(W - W_appr)/norm(W)

% How good is cross approximation
%dif_a = norm(expW_cross - expW)
%dif_r = norm(expW_cross - expW)/norm(expW)
%====== Errors are very small -> works correct



% Now Gregory Series to compute f(w)=log(w) . It is not pointwise.
sz=size(W);



A = Id - W
rank(A)
B = Id + W
rank(B)

%Compute the Hadamard inverse of B

alpha = 1/power(norm(B),2);
alpha = 20*alpha;
V0 = alpha*B;
ER = Id - B.*V0;
err_check = norm(ER)
err_check2 = norm(Id - alpha*B.*B)
invHB = inverse_with_truncation(@inverse_psi, V0, B, Id, trunc_eps)
abs_err = norm(Id - B.*invHB)
rel_err = norm(Id - B.*invHB)/n^d

%One more check
fun_inv = @(x) 1./x
invB_cross = funcrs2(B, fun_inv, 1e-12, B ,18);  % use the cross method to compute f(w)
abs_err = norm(Id - B.*invB_cross)
rel_err = norm(Id - B.*invB_cross)/n^d

err3 = invHB - invB_cross;
err3=round(err3, trunc_eps);
norm(err3)
norm(err3)/norm(invB_cross)


C0 = A.*invHB;
rank(C0)

K=4
Y= iteration_with_truncation(K, C0,  Id, trunc_eps); % compute log(C0) here
Y

%Check ourself , compare with Cross approximation
norm(logW_cross - Y)
norm(logW_cross - Y)/norm(Y)

% X - exp(log(X)) = 0
explog_cross = funcrs2(Y, fun_exp, 1e-12, W ,18);  % use the cross method to compute f(w)
display(explog_cross)
norm(explog_cross - W)
norm(explog_cross - W)/norm(W)



function y = iteration_with_truncation(K, C0,  Id, eps)
  sum = C0;
  for k = 1:K
    prod = Id;
    for j=1:(2*k+1)
        prod = prod .* C0;
        sz=size(prod)
        prod.r
        prod = round(prod, eps);
        prod.r
        sz=size(prod)
    end    
   sum = sum + prod/(2*k+1);
   %rel_err = norm(exact - (-2)*sum)/norm(exact)
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

% fun = @(x) sqrt(x);
% tt=funcrs2(lp1, fun, 1e-12, lp1, 8);
% 

%e2= Id - A.*V

function y  = inverse_with_truncation(function_iteration, xV, xA, I, eps)

  err=10000;
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


function Y = inverse_psi(V, A, I)
   Y = V.*(2*I - A.*V);
end
