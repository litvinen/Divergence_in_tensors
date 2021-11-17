% Compute analytical formula for the Hellinger distance

clc;
clear;
format long

d = 8;
n = 128;
r = 0.5;
tol = 1e-14;
% Spatial domain
a = 20;

h = 2*a/n;
x = ((-n/2+1):(n/2))'*h;
volume = ((max(x) - min(x))/n)^d;
mu1 = 0*ones(1,d);
mu2 = 0*ones(1,d);
sigm1 = 1;
sigm2 = 5;
S1 = sigm1^2*eye(d);
iS1=inv(S1);
S2 = sigm2^2*eye(d);
iS2=inv(S2);
detS1 = det(S1);
detS2 = det(S2);
%rng('default')  % For reproducibility
Temp1 = (r*iS1 + (1-r)*iS2);
%KK = power(det(Temp1), -1/2);

%KK = KK/(power(detS1,r/2) * power(detS2, (1-r)/2))


[U S V] = svd(Temp1);
Temp2 = sum(log(diag(S)));
K = -0.5 * Temp2;
[U S V] = svd(S1);
Temp3 = sum(log(diag(S)))*r/2;
[U S V] = svd(S2);
Temp4 = sum(log(diag(S)))*(1-r)/2;

KKK = K - (Temp3+Temp4);
%KKK = exp(K)
QF = (mu1 - mu2) * inv(r*S2 + (1-r) *S1) * (mu1 - mu2)'
KKK = KKK  + (r*(r-1) * QF / 2)
%KK = KK*exp(r*(r-1) * QF / 2)
K=exp(KKK)

DH2sq  = (1-K) 
DH  = power ((1-K), 0.5) 




g0 = tt_unit(n,1,n/2);
g0 = mtkron(repmat({g0}, 1, d));


fun_P = @(x)exp(-sum((x-mu1(:)').^2,2)/(2*sigm1^2))/ ((2*pi)^(d/2)*sqrt(detS1));
fun_Q = @(x)exp(-sum((x-mu2(:)').^2,2)/(2*sigm2^2))/ ((2*pi)^(d/2)*sqrt(detS2));
fun_sqHellinger = @(x) (sqrt(fun_P(x)) - sqrt(fun_Q(x))).*(sqrt(fun_P(x)) - sqrt(fun_Q(x)));
tic
sqHellinger_tensor = amen_cross_s(n*ones(d,1), @(ind)fun_sqHellinger (reshape(x(ind), [], d)), tol, 'vec',true, 'y0', g0);
display(sqHellinger_tensor)
dt=whos('sqHellinger_tensor'); 
MB=dt.bytes*9.53674e-7
sqHellinger_value = dot(tt_ones(n,d)*h^d, sqHellinger_tensor)
sqHellinger_value = sqHellinger_value/2;
Hellinger_TT = sqrt(sqHellinger_value)
comptime = toc
ea = abs(DH - Hellinger_TT)
er = abs(DH - Hellinger_TT)/DH




