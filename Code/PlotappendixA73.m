clear; clc;
%Program 2: Plotting MAE vs number of samples 
%for interval [0,1]
% function y=exp(-x)
% a=0, b=1, therefore interval is ignored 
exact=0.63212; % exact solution of the function
dim=60; % dimension=60 
rng default 
for n=10000:10000:100000 % samples
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
        xHaltonr=haltr(:, d); % get x vector 
        xSobol=sobol(:, d); % get x vector 
        xSobolr=sobolr(:, d); % get x vector 
        yHalton=exp(-xHalton); % find f(x) 
        yHaltonr=exp(-xHaltonr); % find f(x) 
        ySobol=exp(-xSobol); % find f(x) 
        ySobolr=exp(-xSobolr); % find f(x) 
        yyHalton(d)=mean(yHalton);
        yyHaltonr(d)=mean(yHaltonr);
        yySobol(d)=mean(ySobol);
        yySobolr(d)=mean(ySobolr);
    end
    maeHalton(n/10000)=mean(abs(yyHalton-exact)); % MAE for Halton sequence
    maeHaltonr(n/10000)=mean(abs(yyHaltonr-exact)); % MAE for Halton sequence
    maeSobol(n/10000)=mean(abs(yySobol-exact)); % MAE for Sobol sequence
    maeSobolr(n/10000)=mean(abs(yySobolr-exact)); % MAE for Sobol sequence
end
p(1) = plot(1:10,maeHalton, ':*k', 'linewidth', 2);
hold on
p(1) = plot(1:10,maeHaltonr, ':*g', 'linewidth', 1);
p(3) = plot(1:10, maeSobol, '-r', 'linewidth', 1);
p(4) = plot(1:10, maeSobolr, '--b', 'linewidth', 1);
hold off
axis([1 10 0 0.0004])
legend('Halton','Halton Randomized','Sobol','Sobol Randomized')
grid ON
xlabel('Nx10000')
ylabel('MAE')
title('f(x)=exp(-x)')