% Created on Oct 15, 2021
% Assumed \pcf is given, then \pdf is computed via FFTn + TTtensors
% And then log(\pdf) - the Entripy is computed
clc;
clear;


format long;


d = 16;
n = 16;
Id=tt_ones(n,d)/n^d;
%a = 128;
a=2*pi;
x = (-a+2*a/n:2*a/n:a)';


tic
freq = (-pi/(x(2)-x(1))+2*pi/(2*a):2*pi/(2*a):pi/(x(2)-x(1)))';

fft_magn = n*sqrt(2*pi)/(2*a);
% ifft_magn = fft_magn / n

P = sparse(n,n);
P(1:n/2, n/2:-1:1) = speye(n/2);
P(n/2+1:n, n:-1:n/2+1) = speye(n/2);

tol = 1e-12;
alpha=1.5; % working for alpha={0.7, 1.0, 1.5, 2.0}
mu = 0.0.*(1:d)';
    
%gfun = @(x)exp(-x(:,1).^alpha/2-x(:,2).^alpha/2-x(:,3).^alpha/2) / (2*pi)^(d/2);
%gfun = @(x)exp(-sqrt(x(:,1).^2/2+x(:,2).^2/2+x(:,3).^2/2+x(:,4).^2/2+x(:,5).^2/2).^alpha) / (2*pi)^(d/2);
gfun_P = @(x)exp(1i*x*mu -sqrt(sum(x.^2,2)/2).^alpha) / (2*pi)^(d/2);
alpha = 0.9;
gfun_Q = @(x)exp(1i*x*mu -sqrt(sum(x.^2,2)/2).^alpha) / (2*pi)^(d/2);

g0 = tt_unit(n,1,n/2);
g0 = mtkron(repmat({g0}, 1, d));

g_P = amen_cross_s(n*ones(d,1), @(ind)gfun_P(reshape(x(ind), [], d)), tol, 'vec',true, 'y0', g0);
g_Q = amen_cross_s(n*ones(d,1), @(ind)gfun_Q(reshape(x(ind), [], d)), tol, 'vec',true, 'y0', g0);
    
%    g = amen_cross_s(n*ones(d,1), @(ind)gfun(reshape(x(ind), [], d)), tol, 'vec',true);
    
f_P = g_P;
for i=1:d
        fi = f_P{i}; % dimensions r1 x n x r2
        [r1,~,r2] = size(fi);
        fi = permute(fi, [2,1,3]);
        fi = reshape(fi, n, r1*r2);
        fi = P*fi;
        fi = ifft(fi);
        fi = fi/(fft_magn/n);
        fi = P*fi;
        fi = reshape(fi, n, r1, r2);
        fi = permute(fi, [2,1,3]);
        f_P{i} = fi;
end

f_Q = g_Q;
for i=1:d
        fi_Q = f_Q{i}; % dimensions r1 x n x r2
        [r1,~,r2] = size(fi_Q);
        fi_Q = permute(fi_Q, [2,1,3]);
        fi_Q = reshape(fi_Q, n, r1*r2);
        fi_Q = P*fi_Q;
        fi_Q = ifft(fi_Q);
        fi_Q = fi_Q/(fft_magn/n);
        fi_Q = P*fi_Q;
        fi_Q = reshape(fi_Q, n, r1, r2);
        fi_Q = permute(fi_Q, [2,1,3]);
        f_Q{i} = fi_Q;
end
%ind=(f>0.9);
f_P = f_P*(2*pi/(x(2)-x(1)))^d;
display(f_P)
%xP = real(full(f_P));

f_Q = f_Q*(2*pi/(x(2)-x(1)))^d;
display(f_Q)
%xQ = real(full(f_Q))

dZ_P = dot(Id, f_P) - 1
dZ_Q = dot(Id, f_Q) - 1
if d==1
  plot(x,real(f),'.')
end    
if d==2 
   % mesh(x,x,real(f))
    figure
    surf(x,x,reshape(full(real(f_P)),n,n))
    figure
    surf(x,x,reshape(full(real(f_Q)),n,n))
    
    %surfc(x,x,reshape(full(real(f)),n,n))
    %mesh(full(real(f), [x x]))
    %mesh(full(real(f), [n n]))
end    
%f.*(log(f) - log(g))

% fun_log = @(x) log(x);
% logP = amen_cross_s({f_P}, fun_log, 1e-12);
% logQ = amen_cross_s({f_Q}, fun_log, 1e-12);
% KLD = f_P.*(logP - logQ);



fun_Hellinger = @(x) power((sqrt(x(:,1)) - sqrt(x(:,2))),2);
Hellinger_cross = amen_cross_s({f_P, f_Q}, fun_Hellinger, tol);
display(Hellinger_cross)
dt=whos('Hellinger_cross'); 
MB=dt.bytes*9.53674e-7

% Compute Hellinger distance
Hellinger_value = real(dot(Id, Hellinger_cross))/2;% *(2*pi/(x(2)-x(1)))^d 
Hellinger_value = sqrt(Hellinger_value)
mytime = toc

%in_use = monitor_memory_whos
%norm(KLD - KLD_cross)
%norm(KLD - KLD_cross)/norm(KLD)



%xlogx_cross2 = amen_cross_s({f}, fun_xlog, 1e-12);
%display(xlogx_cross2)

