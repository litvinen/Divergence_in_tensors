% % Example to compute the Kullback-Leibler divergence (KLD) of two Gaussians, 
% % the first mean is 1 and the second -1, both variances are equal to 0.2
% clc;
% clear;
% h=0.1; % discretisation step size
% % Define Interval [-10,10]
% a=-10;
% b=10; 
% x=[a:h:b]
% y1=normpdf(x,-1,1)
% %plot(x,y1,'*')
% y2=normpdf(x,1,1)
% 
% KLD_1d=h*sum(y1.* (log(y1) - log(y2)) )

%=====
clc;
clear;
d=2;
% 
% N=100;
% mu = [1 -1]; 
% Sigma = [.9 .4; .4 .3];
% [X1,X2] = meshgrid(linspace(-1,3,N)', linspace(-3,1,N)');
% %plot(X1,X2,'.')
% X = [X1(:) X2(:)];
% p = mvnpdf(X, mu, Sigma);
% 
% surf(X1,X2,reshape(p, N, N));

% mu = [1.4 -2.8 5.4]; 
% Sigma = [.9 .4 .5; .4 .3 .34;  .48 .34 .24];
% Sigma = Sigma *Sigma';
% r = mvnrnd(mu, Sigma, 15000);
% for i=1:3
%     sum(r(:,i))/15000
% end    
% plot3(r(:,1),r(:,2),r(:,3),'.');

N=10;
M=10; %randomly sample M points
mu1 = 0*ones(1,d);
mu2 = 0*ones(1,d);
Sigma1 = eye(d);
Sigma2 = 10*eye(d);
rng('default')  % For reproducibility

KLD_anal_12 = 0.5*(trace(inv(Sigma2)*Sigma1) + (mu2-mu1)*inv(Sigma2)*(mu2-mu1)' - d + log(det(Sigma2)/det(Sigma1)))
KLD_anal_21 = 0.5*(trace(inv(Sigma1)*Sigma2) + (mu1-mu2)*inv(Sigma1)*(mu1-mu2)' - d + log(det(Sigma1)/det(Sigma2)))

%X1 = mvnrnd(mu1,Sigma1,M);
%X2=X1;
%X2 = mvnrnd(mu2,Sigma2,M);
%Y1=mvnpdf(X1, mu1, Sigma1);
%Y2=mvnpdf(X2, mu2, Sigma2);
%scatter3(X1(:,1),X1(:,2),Y1,'.')
%hold all
%scatter3(X2(:,1),X2(:,2),Y2,'.')
%xlabel('X1')
%ylabel('X2')
%zlabel('Probability Density')
% 
%Sample1 = mvnrnd(mu1, Sigma1, M);
%Sample2 = mvnrnd(mu2, Sigma2, M);
h=0.1;
spaceX= -10:h:10;
spaceY= -10:h:10;
nx = size(spaceX,2);
ny = size(spaceY,2);
onesX=ones(ny,1);
onesY=ones(nx,1);
XX1=kron(spaceX, onesY);
XX1=reshape(XX1,nx,ny);

YY1=kron(onesY, spaceY');
YY1=reshape(YY1,nx,ny);

N=size(spaceX,2);
[X1,X2] = meshgrid(spaceX, spaceY');
X = [X1(:) X2(:)];
normX = X1-XX1;
normX = norm(X1-XX1)
normY = X2-YY1;
normY = norm(X2-YY1)

XX = [XX1(:) YY1(:)];
norm(X-XX)
% plot(X1(:), X2(:), '.')
Y1 = mvnpdf(XX, mu1, Sigma1);
Y2 = mvnpdf(XX, mu2, Sigma2);
%scatter3(X1(:), X2(:), Y1, '.')
%hold all
%Y2 = mvnpdf(Sample2, mu2, Sigma2);
%scatter3(X1,X2, (Y2),'.')



logY1 = log(Y1);
logY2 = log(Y2);
dif = Y1.*(logY1 - logY2);
KLD_num = sum(dif)*power(h,2)




%S2=mvnrnd(mu2, Sigma2, M);
%figure
%plot(S1)
%hold all
%plot(S2)
% for i=1:M
%   scatter3(X1(:,1),X1(:,2),S1(i),'.')
%   hold all
%   scatter3(X2(:,1),X2(:,2),S2(i),'.')
%   hold all
% end  
%KLD=sum(S2.* (log(S2) - log(S1)) )
%KLD=sum(S1.* (log(S1) - log(S2)) )
