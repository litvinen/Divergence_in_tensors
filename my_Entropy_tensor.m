% Created on Oct 15, 2021
% Assumed \pcf is given, then \pdf is computed via FFTn + TTtensors
% And then log(\pdf) - the Entripy is computed
clc;
clear;



format long;

d = 5;
n = 16;

%a = 128;
a=2*pi;
x = (-a+2*a/n:2*a/n:a)';



freq = (-pi/(x(2)-x(1))+2*pi/(2*a):2*pi/(2*a):pi/(x(2)-x(1)))';

fft_magn = n*sqrt(2*pi)/(2*a);
% ifft_magn = fft_magn / n

P = sparse(n,n);
P(1:n/2, n/2:-1:1) = speye(n/2);
P(n/2+1:n, n:-1:n/2+1) = speye(n/2);

tol = 1e-6;
alpha=0.5; % working for alpha={0.7, 1.0, 1.5, 2.0}
mu = 0.0.*(1:d)'
    
%gfun = @(x)exp(-x(:,1).^alpha/2-x(:,2).^alpha/2-x(:,3).^alpha/2) / (2*pi)^(d/2);
%gfun = @(x)exp(-sqrt(x(:,1).^2/2+x(:,2).^2/2+x(:,3).^2/2+x(:,4).^2/2+x(:,5).^2/2).^alpha) / (2*pi)^(d/2);
gfun = @(x)exp(1i*x*mu -sqrt(sum(x.^2,2)/2).^alpha) / (2*pi)^(d/2);

g0 = tt_unit(n,1,n/2);
%t_=full(g0)
g0 = mtkron(repmat({g0}, 1, d));
%myfun = @(ind)gfun(reshape(x(ind), [], d));
g = amen_cross_s(n*ones(d,1), @(ind)gfun(reshape(x(ind), [], d)), tol, 'vec',true, 'y0', g0);
    
%    g = amen_cross_s(n*ones(d,1), @(ind)gfun(reshape(x(ind), [], d)), tol, 'vec',true);
    
f = g;
for i=1:d
        fi = f{i}; % dimensions r1 x n x r2
        [r1,~,r2] = size(fi);
        fi = permute(fi, [2,1,3]);
        fi = reshape(fi, n, r1*r2);
        fi = P*fi;
        fi = ifft(fi);
        fi = fi/(fft_magn/n);
        fi = P*fi;
        fi = reshape(fi, n, r1, r2);
        fi = permute(fi, [2,1,3]);
        f{i} = fi;
end

%ind=(f>0.9);

f= real(f*(2*pi/(x(2)-x(1)))^d);
dZ = dot(tt_ones(n,d)/n^d, f) - 1
if d==1
  plot(x,real(f),'.')
end    
if d==2 
   % mesh(x,x,real(f))
    surf(x,x,reshape(full(real(f)),n,n))
    
    %surfc(x,x,reshape(full(real(f)),n,n))
    %mesh(full(real(f), [x x]))
    %mesh(full(real(f), [n n]))
end    


fun_sqrt = @(x) sqrt(x);
sqrt_cross = funcrs2(f, fun_sqrt, 1e-12, f ,18);  % use the cross method to compute f(w)
display(sqrt_cross)
sqrt_cross  = round(sqrt_cross, 1e-2 )
display(sqrt_cross)
norm(f - sqrt_cross.*sqrt_cross)/n^d


fun_log = @(x) log(x);
logx_cross = funcrs2(f, fun_log, 1e-12, f ,18);  % use the cross method to compute f(w)
display(logx_cross)

fun_xlog = @(x) x.*log(x);
xlogx_cross = funcrs2(f, fun_xlog, 1e-12, f ,18);  % use the cross method to compute f(w)
display(xlogx_cross)

logx_cross2 = amen_cross_s({f}, fun_log, 1e-12);
display(logx_cross2)

xlogx_cross2 = amen_cross_s({f}, fun_xlog, 1e-12);
display(xlogx_cross2)

% Comparison of two methods
norm(logx_cross2 - logx_cross)/norm(logx_cross)
norm(xlogx_cross2 - xlogx_cross)/norm(xlogx_cross)
%

% Compute Entropy
Id = tt_ones(n,d)/n^d;
entropy = real(dot(Id, xlogx_cross))% *(2*pi/(x(2)-x(1)))^d 
