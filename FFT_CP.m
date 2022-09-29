%
% My implementation of fast fourier transform in the Heston model
%
% Author: Raza Ali Mahmood
% Date: 21th April 2021
% 
%
% FFT Option Pricing, we have taken our inspiration from 
% Option Valuation Using the Fast Fourier Transform, 
%       Peter Carr, March 1999, pp 10-11
%
%
% Input:   model    - Takes in the name of the model we are using ti the
%                       appripiate characteristic function
%          S_0      - Spot price of the underlying asset
%          K        - Strike price
%          T        - Maturity, time measured in (years)
%          r        - The risk-free rate
%          d        - Dividend yield 

%          varagin  - s an input variable in a function definition
%                     statement that enables the function to accept any
%                     number of input arguments. Specify varargin by using lowercase
%                     characters. After any explicitly declared inputs, include varargin 
%                     as the last input argument 


% Output:  Price - The price of the European call option.
%

function price = FFT_CP(model,n,S_0,K,T,r,d,varargin)
% Defining basic variable for later use.
logS = log(S_0);
logK = log(K);

%option_Alpha= optimalAlpha(model,lnS,lnK,T,r,d,varargin{:});
option_Alpha = .75;
Discount_Factor = exp(-r*T);


% Predefined parameters
PFFT_N = 2^n;                               % Must be a power of two (2^14)

% Spacing of psi integrand
FFT_beta = 0.05;                            

% Spacing for log strike output equation (23)
FFT_lambda = (2 * pi) / (PFFT_N * FFT_beta);  

% Follows from equation (20)
FFT_b = (PFFT_N * FFT_lambda) / 2;           

% prepping for equation 19
uvec = 1:PFFT_N;
% Equation (20), gives us the log strike levels ranging from -b to b where:
nu = - FFT_b + FFT_lambda * (uvec - 1);   

% Loop for FFT step
jvec = 1:PFFT_N;
vj = (jvec-1) * FFT_beta;

% Applying FFT
temporary = Discount_Factor * psi(model,vj,option_Alpha,logS,T,r,d,varargin{:}) .* exp(1i * vj * (FFT_b)) * FFT_beta;
% Applying simpson's rule
temporary = (temporary / 3) .* (3 + (-1).^jvec - ((jvec - 1) == 0) );  

% Call price vector resulting in equation (24)
call_price_vector = real(exp(-option_Alpha .* nu) .* fft(temporary) / pi);        

indexOfStrike = floor((logK + FFT_b)/FFT_lambda + 1); 
iset = max(indexOfStrike)+1:-1:min(indexOfStrike)-1;
xp = nu(iset);
yp = call_price_vector(iset);
price = real(interp1(xp,yp,logK));
end

% Analytical formula for psi in equation ( 6 ) of Madan's paper
% Here we use the Characteristic library inspired from
% Kienitz, Joerg and Wetterau, Daniel to make the code nicer
function ret = psi(model,v,alpha,varargin)
  ret = exp(feval(@CFL, model, v - (alpha + 1) * 1i,varargin{:})) ./ (alpha.^2 + alpha - v.^2 + 1i * (2 * alpha + 1) .* v);
end
