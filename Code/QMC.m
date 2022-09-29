%
% My implementation of quasi-Monte Carlo Simulation and randomized
% quasi-Monte Carlo with broadie kaya shceme in the Heston model
% to price European call option 
% 
% Author: Raza Ali Mahmood
% Date: 21th April 2021
% 
% Prepping is written when pre-allocating size of vectors for variables
% This allows for less error and faster computation
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
%          Time  - The time estimate of how long the Monte Carlo
%                  simulation took till finish
%

function [Price, Time, pathS] = QMC(S0,r,T,K,V_0,theta,kappa,sigma,rho,Nsimp,CallorPut)
% Setting up my time estimator
Timer = clock;
NTime = 1;
% Sobol sequence setup, with randomization, skip and leap 
qsobolset = sobolset(1,'Skip',100,'Leap',10);   
%qsobolset = scramble(qsobolset,'MatousekAffineOwen');

% Use either pseudorandom sequence rand or quasirandom sequence sobol
%LD = rand(Nsimp,NTime);              % use for random numbers
LD = net(qsobolset,Nsimp);            % use for sobol numbers

% Creating my stepsize and the sizing of my grid parameters
dt = T;                                              
N = 800; % choose to cut it into 3200 intervals from 0 to 200
int = 0.25;      
grid = int * (1:N);   

% Creating tolerance level
options = optimset('TolX',1e-8);                     
                                         

% Prepping size of path Variance, Integrated Variance and defining variance
pathV = zeros(Nsimp,NTime+1); pathV(:,1) = V_0;     
Vold = V_0*ones(Nsimp,1);
I_var = zeros(Nsimp,1);                         
Variance = sigma^2;

% Prepping size of path S
pathS = zeros(Nsimp,NTime+1); 
pathS(:,1) = S0;        


% Analytical properties C_0 psi and lambda 
psi = 4 * kappa * theta / Variance;
lambda = 4 * kappa * exp(-kappa*dt) / (Variance * (1 - exp(-kappa*dt)));
C_0 = Variance * (1 - exp(-kappa*dt)) / (4 * kappa);

% Defining variable for the characteristic function, A(a), B, C(a),
% Gamma(a). We defined more variables for more ease
gammaa = sqrt(kappa^2 - 2 * Variance * 1i * grid );

% Defining phi construct
phi1 = gammaa .* exp( -0.5 * (gammaa-kappa) * dt ) ...
    *(1 - exp(-kappa*dt)) ./ (kappa * (1 - exp(-gammaa*dt)));

% Defining some of C
C2 = kappa * (1 + exp(-kappa * dt)) / (1 - exp(-kappa * dt));
C3 = gammaa .* (1 + exp(-gammaa * dt)) ./ (1-exp(-gammaa * dt));


% Constants used in A, B
C0const = 2 ./ grid  / pi * int;
A1const = 4 .* gammaa .* exp(-0.5 * gammaa * dt) ...
    ./ ( Variance * (1 - exp(-gammaa * dt)) );

B1construct = 4 * kappa * exp(-0.5 * kappa * dt) ...
    / (Variance * (1 - exp(-kappa * dt)));


% We can optimize the speed of computation by calculating the difference in
% C(a) outside the loop, as they are two constants C2 and C3 and defining
% the value for the bessel function indicatore psi/2-1 as mu
CDiff = C2-C3;         
mu = 0.5 * psi - 1;       

for col = 2:NTime+1    
    %First step of the algorithm is calculating V_t (Vnew) conditioning on
    % V_s (V_old)
    V_s = lambda * Vold;
    Vnew = C_0 * ncx2rnd(psi,V_s,Nsimp,1); % Using the built in non-central
    % chi-squared distribution
    
    % Prepping for the second step of the algorithm
    % Calculating the variables on page 4 A, B, and C
    C1 = (Vold + Vnew) / (Variance);   
    B1 = besseli(mu,sqrt(Vold .* Vnew) * B1construct);
    A1repmat = besseli(mu,sqrt(Vold(:) .* Vnew(:)) * A1const);  
    
    % Initizilization value for the grid search function fzero
    % It is not described in the algorithm, but many cases of root
    % searching it is nessecary to give a starting point 
    % Simply chose the discounted value of the V_s
    startval = theta + (Vold - theta) * exp(-kappa*dt); % start for fzero
     
    % Loop for the second step
    for row = 1:Nsimp
        % Calculating the charateristic function Phi
        CharFunc = real(phi1 .* exp(C1(row) .* CDiff) ...
            .* A1repmat(row,:) ./ B1(row));             
        
        % Begin the grid search on the G(x), note that we approximate the
        % integral with a sum using our quasirandom numbers as a uniform
        % variable using the starting value to increase the speed of the 
        % grid search functino
        I_var(row) = fzero(@(x) int*x/pi ...
            + sum(C0const .* sin(grid * x) .* CharFunc)  ...
          - LD(row,col-1),startval(row),options);      
        % Using what we know from the full truncation scheme to fix 
        % negative values in our variance
        I_var(row) = max(I_var(row),0);       
    end
    
    % Step 3, fairly simple use equation (2.17) was what it was last time
    % Sampled integrated volatility
    SI_Vol = (1/sigma)*(Vnew - Vold - kappa * theta * dt + kappa * I_var);    
    
    % Generate an independent sample from a standard normal distribution
    % using the integrated variance and sampled integrated vol
    mu = log(pathS(:,col-1)) + r * dt - 0.5 * I_var + rho * SI_Vol;  
    
    % prepping for equation (2.18)
    sig = (1 - rho^2) * I_var;            
    
    % Step 4 creating the path for the stock S, using the exact solution
    pathS(:,col) = exp(mu + randn(Nsimp,1) .* sqrt(sig));
    
    % prepping for the next iteration with the new updated V becoming the
    % old V, and also saving the paths of V

    Vold = Vnew; 
    pathV(:,col) = Vnew;
end
% Last step calculate the price of the paths of S, here we used our built
% home made Call_put option, which simply takes the mean of the max(S-K,0)
% And at last the code is dont and we save the time estimate.
Price = Call_Put(pathS,K,CallorPut);
Time = etime(clock,Timer);