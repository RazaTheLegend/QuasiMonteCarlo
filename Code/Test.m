clear; clc;
% Parameters
S_0 = 100;    % Spot prices
K = 100;    % Strike price
T = 1;      % maturities
r = 0;      % discount factors
d = 0;      % dividends
% Specify the model 
V_0 = 0.04;      % instantanuous variance of base parameter set  
theta = 0.04;      % long term variance of base parameter set
kappa = 0.25;      % mean reversion speed of variance of base parameter set
sigma = 0.5;       % volatility of variance of base parameter set
rho = -0.6;        % correlation of base parameter set

% Simulation parameters
Nsimp = 100000; CallorPut=1; stepsize=100;
model = 'Heston';
n=19;

Price =ones(5,2);
% method 1 FFT using myheston
Price(1,2)=FFT_CP(model,n,S_0,K,T,r,d,V_0, theta, kappa, sigma, rho);

% method 1 EMC using myheston
Price(2,2)=MyHeston(S_0,r,T,K,V_0,theta,kappa,sigma,rho,Nsimp,stepsize);

% method 2 BKMC with 
Price(3,2)=BKMC(S_0,r,T,K,V_0,theta,kappa,sigma,rho,Nsimp,CallorPut);

% method 3  QMC with 
Price(4,2)=QMC(S_0,r,T,K,V_0,theta,kappa,sigma,rho,Nsimp,CallorPut);

% method 4 RQMC with 
Price(5,2)=RQMC(S_0,r,T,K,V_0,theta,kappa,sigma,rho,Nsimp,CallorPut);
 



