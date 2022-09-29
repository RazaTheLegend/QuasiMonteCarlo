clear; clc;
% Code need for the cascade of plots comparing the speed of convergence and
% different sizes of sample path and speed estimates.
% Simple method to observe the progess in a loop or similar fashion
f = waitbar(0, 'Starting');
% Parameters
S_0 = 100;    % Spot prices
T = 1;      % Time till maturity 1, meaning 1 year.
K = 100;    % Strike price
r = 0;      % discount factors
d = 0;      % dividends
% Specify the model 
V_0 = 0.04;      % instantanuous variance of base parameter set  
theta = 0.04;      % long term variance of base parameter set
kappa = 0.25;      % mean reversion speed of variance of base parameter set
sigma = 0.5;       % volatility of variance of base parameter set
rho = -0.6;        % correlation of base parameter set


% Simulation parameters
CallorPut = 1; % Call = 1 -> Call; Call = 0 -> Put

NSim = 100000;
Price1 = ones(1,NSim);
Time1 = ones(1,NSim);
Price2 = ones(1,NSim);
Time2 = ones(1,NSim);
Price3 = ones(1,NSim);
Time3 = ones(1,NSim);
Price4 = ones(1,NSim);
Time4 = ones(1,NSim);
Price11 = ones(1,NSim);
Price21 = ones(1,NSim);
Price31 = ones(1,NSim);
spacing = linspace(1,NSim,1);
x = CallPricingFFT('Heston',18,S,K,T,r,d,v, theta, kappa, omega, rho);
   
for i = spacing
rng(i)  
% method 1 creating paths and passing to payoff
[Price1(i), Time1(i)] = RQMC(S_0, r, T, K, V_0, theta, kappa, sigma, rho ,i, CallorPut);
[Price2(i), Time2(i)] =  QMC(S_0, r, T, K, V_0, theta, kappa, sigma, rho ,i, CallorPut);
[Price3(i), Time3(i)] = BKMC(S_0, r, T, K, V_0, theta, kappa, sigma, rho ,i, CallorPut);
[Price4(i), C_F, Time4(i)] = MyHeston(S_0, r, T, K, V_0, theta, kappa, sigma, rho, 100, i,4);

    waitbar(i/NSim, f, sprintf('Progress: %d %%', floor(i/NSim*100)));
    pause(0.1);
end
close(f)
X = spacing;

%Plot for option price
figure;
plot(X,Price1(spacing),'m-',X,Price2(spacing),'g:', ...
    X,price3(spacing),'b--',X,price4(spacing),'r--')
yline(6.3219,'-','True Price')
legend('RQMC','QMC','BKMC','EMC','')
title("Convergence")
xlabel("Number of Simulations")
ylabel("Price")
axis([0 NSim 6.1 6.5])


%Plot for time
figure;
plot(X,Time1(spacing),'m-',X,Time2(spacing),'g:', ...
    X,Time3(spacing),'b--',X,Time4(spacing),'r--')
legend('RQMC','QMC','BKMC','EMC')
title("Speed of Simulation")
xlabel("Number of Simulations")
label("Time")