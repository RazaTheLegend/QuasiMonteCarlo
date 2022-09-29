clear; clc;
%Program 1: Plotting absolute error vs dimensions
%for interval=[0,1]
% function y=sin(2x)
% a=0, b=1, therefore interval is ignored 
exact=0.70807;
dim=50; % for 60 dimensions
n=100000; % number of samples
rng default
h=haltonset(dim);
halt=net(h, n);
hr=scramble(h,'RR2');
haltr=net(hr, n);
s=sobolset(dim);
sobol=net(s, n);
sr=scramble(s,'MatousekAffineOwen');
sobolr=net(sr, n);
for d=1:dim
    xHalton=halt(:, d); % get x vector 
    yHalton=sin(2*xHalton); % find f(x)

    xHaltonr=haltr(:, d); % get x vector 
    yHaltonr=sin(2*xHaltonr); % find f(x)

        xSobol=sobol(:, d);% get x vector
        ySobol=sin(2*xSobol); % find f(x)

        xSobolr=sobolr(:, d);% get x vector
        ySobolr=sin(2*xSobolr); % find f(x)

        yyHalton(d)=mean(yHalton);
        yyHaltonr(d)=mean(yHaltonr);
        yySobol(d)=mean(ySobol); 
        yySobolr(d)=mean(ySobolr);
end
errorHalton=abs(yyHalton-exact); % absolute
% error for Halton sequence
errorHaltonr=abs(yyHaltonr-exact); % absolute
% error for Halton sequence
errorSobol=abs(yySobol-exact); % absolute
% error for Sobol sequence
errorSobolr=abs(yySobolr-exact); % absolute
 % error for Sobol sequence
plot(1:dim, errorHalton, ':.k','linewidth',1)
hold all
plot(1:dim, errorHaltonr, ':.g','linewidth',1)
plot(1:dim, errorSobol, '-*r','linewidth',1,'MarkerSize',2.5)
plot(1:dim, errorSobolr, '-ob','linewidth',1,'MarkerSize',2)
axis([0 dim 0 0.0001])
grid ON
legend('Halton','Halton Randomized','Sobol','Sobol Randomized')
xlabel('Dimension')
legend('Halton','Halton Randomized','Sobol','Sobol Randomized')
grid ON
ylabel('Absolute Error')
title('f(x)=sin(2x)')
