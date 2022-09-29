%
% My implementation of Monte Carlo Simulation with Euler discretization
% with Full Truncation in the Heston model to price European call option
% 
% Author: Raza Ali Mahmood
% Date: 21th April 2021
% 
% Input: S_0      - Spot price of the underlying asset
%        K        - Strike price
%        T        - Maturity, time measured in (years)
%        r        - The risk-free rate
%        V_0      - Initial variance
%        theta    - Long-term average variance
%        kappa    - Mean reversion speed of variance of base parameter se
%        sigma    - Volatility of variance
%        rho      - Instantaneous correlation coefficient
%        stepsize - Number of steps per year
%        Nsimp    - Number of simulation paths

% Output:  Price - The price of the European call option.
%          C_F   - The cashflow from each sample path
%          Time  - The time estimate of how long the Monte Carlo
%                  simulation took till finish
%
function [Price, C_F, Time] = MyHeston(S_0,r,T,K,V_0,theta,kappa,sigma,rho,Nsimp,stepsize)

f_1 = @(x)x;
f_2 = @(x)max(0,x); 
f_3 = @(x)max(0,x);

dt = 1/stepsize; % Step per year will give us the stepsize dt

% Initialize matrices to record values during the simulation
V_tilde = zeros(Nsimp , T/dt+1);
V       = zeros(Nsimp , T/dt+1);
ln_S_c  = zeros(Nsimp , T/dt+1); % for S0

% Set values at t=0
V_tilde(:,1) = V_0;
V(:,1)       = V_0;
ln_S_c(:,1)    = log(S_0);

% To Generate values of the correlated Brownian motions
dW_V = normrnd(0,sqrt(dt),Nsimp ,T/dt);
dZ   = normrnd(0,1,Nsimp ,T/dt);
dW_S = rho.*dW_V + sqrt(1-rho^2).*dZ.*sqrt(dt);

% The Monte Carlo Simulation loop begins
tic;
for j = 1:T/dt
    V_tilde(:,j+1) = f_1(V_tilde(:,j)) ...
                      - kappa.*dt.*(f_2(V_tilde(:,j))-theta) ...
                      + sigma.*sqrt(f_3(V_tilde(:,j))).*dW_V(:,j);
    V(:,j+1)       = f_3(V_tilde(:,j+1));
    ln_S_c(:,j+1)  = ln_S_c(:,j) + (r-0.5.*V(:,j)).*dt ...
                      + sqrt(V(:,j)).*dW_S(:,j);

end
Time=toc;
% Revert log(S) to S.
S_T_c = exp(ln_S_c(:,end));

% Calculate the option's price
Price = exp(-r*T)*mean(max(S_T_c-K,0));
C_F = exp(-r*T)*(max(S_T_c-K,0));

end