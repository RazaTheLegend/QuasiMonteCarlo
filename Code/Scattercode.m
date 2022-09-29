% % % % % % % % % % % % % % 
% random scatter plot code
% % % % % % % % % % % % % % 
rng default
X = rand(1000,2);
figure;
scatter(X(:,1),X(:,2),4,'red','filled')
axis square;
xlabel(dx);
ylabel(dy);
title('{\bf Uniform Random Scatter}')
% % % % % % % % % % % % % % 
% % % % % % % % % % % % % % 
% sobol scatter plot code
% % % % % % % % % % % % % % 
% % % % % % % % % % % % % % 
rng default
dx=1;
dy=2;
sset=sobolset(100);
x0=net(sset,1000);
size(x0)
figure;
scatter(x0(:,dx),x0(:,dy),4,'b','filled')
axis square;
xlabel(dx);
ylabel(dy);
title('{\bf Sobol Scatter}')

% high dimension sobol scatter plot code
rng default
dx=69;
dy=70;
sset=sobolset(100);
x0=net(sset,4096);
size(x0)
figure;
scatter(x0(:,dx),x0(:,dy),4,'b','filled')
axis square;
xlabel(dx);
ylabel(dy);
title('{\bf High Dimensional Sobol Scatter}')

% high dimension owen scrambled sobol scatter plot code
rng default
dx=1;
dy=2;
sset=sobolset(100);
sset=scramble(sset,'MatousekAffineOwen');
x0=net(sset,1000);
size(x0)
figure;
scatter(x0(:,dx),x0(:,dy),4,'b','filled')
axis square;
xlabel(dx);
ylabel(dy);
title('{\bf Scrambled Sobol Scatter}')
% % % % % % % % % % % % % % 
% % % % % % % % % % % % % % 
% % % % % % % % % % % % % % 
% % % % % % % % % % % % % % 

% high dimension Matousek scrambled skip sobol scatter plot code
rng default
dx=69;
dy=70;
sset=sobolset(100,'skip',10);
sset=scramble(sset,'MatousekAffineOwen');
x0=net(sset,4096);
size(x0)
figure;
scatter(x0(:,dx),x0(:,dy),4,'b','filled')
axis square;
xlabel(dx);
ylabel(dy);
title('{\bf  Scrambled and Skip (10) Sobol Scatter}')
 
rng default
dx=69;
dy=70;
sset=sobolset(100,'skip',10000);
sset=scramble(sset,'MatousekAffineOwen');
x0=net(sset,4096);
size(x0)
figure;
scatter(x0(:,dx),x0(:,dy),4,'b','filled')
axis square;
xlabel(dx);
ylabel(dy);
title('{\bf  Scrambled and Skip (10000) Sobol Scatter}')


% high dimension Matousek scrambled and leap sobol scatter plot code
rng default
dx=69;
dy=70;
sset=sobolset(100,'Leap',10);
sset=scramble(sset,'MatousekAffineOwen');
x0=net(sset,4096);
size(x0)
figure;
scatter(x0(:,dx),x0(:,dy),4,'b','filled')
axis square;
xlabel(dx);
ylabel(dy);
title('{\bf Scrambled and Leap (10) Sobol Scatter}')

rng default
dx=69;
dy=70;
sset=sobolset(100,'Leap',1e2);
sset=scramble(sset,'MatousekAffineOwen');
x0=net(sset,4096);
size(x0)
figure;
scatter(x0(:,dx),x0(:,dy),4,'b','filled')
axis square;
xlabel(dx);
ylabel(dy);
title('{\bf Scrambled and Leap (1e2) Sobol Scatter}')
% % % % % % % % % % % % % % % % % % % % % % % % 
% % % % % % % % % % % % % % % % % % % % % % % % 
% % % % % % % % % % % % % % % % % % % % % % % % 

% high dimension Matousek scrambled leap sobol scatter plot code

rng default
dx=69;
dy=70;
sset=sobolset(100,'Leap',10,'Skip',100);
sset=scramble(sset,'MatousekAffineOwen');
x0=net(sset,4096);
size(x0)
figure;
scatter(x0(:,dx),x0(:,dy),4,'b','filled')
axis square;
xlabel(dx);
ylabel(dy);
title('{\bf Scrambled, Leap (10) and Skip (100) Sobol sequence}')


% Quasirandom sequence cant in general pass the uniform test (k-test).
% which can be seen below, where 5% of the p-values should fall below 0.05
% In a statistical sense, quasi-random numbers are too uniform to
% pass traditional tests of randomness 
sset = sobolset(1,'Skip',1e3,'Leap',1e2);
sset = scramble(sset,'MatousekAffineOwen');

nTests = 1e5;
sampSize = 50;
PVALS = zeros(nTests,1);
for test = 1:nTests
    x = sset(test:test+(sampSize-1),:);
    [h,pval] = kstest(x,[x,x]);
    PVALS(test) = pval;
end
histogram(PVALS,100)
xlabel('{\it p}-values')
ylabel('Number of Tests')
title('{\bf Test of Randomness for Sobol}')

% Quasirandom sequence cant in general pass the uniform test (k-test).
% which can be seen below, where 5% of the p-values should fall below 0.05
% This is not the case for these 

nTests = 1e5;
sampSize = 50;
PVALS = zeros(nTests,1);
for test = 1:nTests
    x = rand(sampSize,1);
    [h,pval] = kstest(x,[x,x]);
    PVALS(test) = pval;
end

histogram(PVALS,100)
h = findobj(gca,'Type','patch');
xlabel('{\it p}-values')
ylabel('Number of Tests')
title('{\bf Scrambled, Leap (10) and Skip (100) Sobol sequence}')
title('{\bf Test of Randomness for MATLAB function}')