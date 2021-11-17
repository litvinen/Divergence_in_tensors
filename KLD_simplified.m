clc;
clear;
% Inverse FFT f(x) = 1/(2*pi) * \int g(omega) exp(j*omega*x) domega
d = 16;
n = 128;
format long

% Spatial domain
a = 10;
h = 2*a/n;
x = ((-n/2+1):(n/2))'*h;



tol = 1e-2;


% Construct image fun
g0 = tt_unit(n,1,n/2);
g0 = mtkron(repmat({g0}, 1, d));






mu1 = 1.1*ones(1,d);
mu2 = 1.4*ones(1,d);
sigm1 = 1.5;
sigm2 = 22.1;
Sigma1 = sigm1^2*eye(d);
invSigma1=inv(Sigma1);
Sigma2 = sigm2^2*eye(d);
invSigma2=inv(Sigma2);
detS1 = det(Sigma1);
detS2 = det(Sigma2);
rng('default')  % For reproducibility

KLD_analytic_12 = 0.5*(trace(invSigma2*Sigma1) + (mu2-mu1)*invSigma2*(mu2-mu1)' - d + log(detS2/detS1))
KLD_analytic_21 = 0.5*(trace(invSigma1*Sigma2) + (mu1-mu2)*invSigma1*(mu1-mu2)' - d + log(detS1/detS2))


fun_P = @(x)exp(-sum((x-mu1(:)').^2,2)/(2*sigm1^2))/ ((2*pi)^(d/2)*sqrt(detS1));
fun_Q = @(x)exp(-sum((x-mu2(:)').^2,2)/(2*sigm2^2))/ ((2*pi)^(d/2)*sqrt(detS2));
fun_KLD = @(x) fun_P(x).*(log(fun_P(x)) - log(fun_Q(x)));
fun_logP = @(x) log(fun_P(x)) ;
fun_logQ = @(x) log(fun_Q(x)) ;
tic
KLD_tensor = amen_cross_s(n*ones(d,1), @(ind)fun_KLD(reshape(x(ind), [], d)), tol, 'vec',true, 'y0', g0);
display(KLD_tensor)
dt=whos('KLD_tensor'); 
MB=dt.bytes*9.53674e-7
KLD_value = dot(tt_ones(n,d)*h^d, KLD_tensor)
comptime = toc
ea = abs(KLD_analytic_12 - KLD_value)
er = abs(KLD_analytic_12 - KLD_value)/KLD_analytic_12

P = amen_cross_s(n*ones(d,1), @(ind)fun_P(reshape(x(ind), [], d)), tol, 'vec',true, 'y0', g0);
Q = amen_cross_s(n*ones(d,1), @(ind)fun_Q(reshape(x(ind), [], d)), tol, 'vec',true, 'y0', g0);
logP = amen_cross_s(n*ones(d,1), @(ind)fun_logP(reshape(x(ind), [], d)), tol, 'vec',true, 'y0', g0);
logQ = amen_cross_s(n*ones(d,1), @(ind)fun_logQ(reshape(x(ind), [], d)), tol, 'vec',true, 'y0', g0);

dt=whos('P'); 
MB1=dt.bytes*9.53674e-7
dt=whos('logP'); 
MB2=dt.bytes*9.53674e-7
2*(MB1+MB2)

display(P)
display(logP)
temp1 = dot(P, logP)
temp2 = dot(P, logQ)
KLDs = (temp1 - temp2)
KLDs * (20^d)/n^d
KLD_analytic_12 

KLD_analytic_12 / KLDs
KLDs / KLD_analytic_12 

