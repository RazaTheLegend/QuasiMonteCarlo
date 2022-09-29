%
%Characteristic Function Library o
%

function phi = CFL(model,u,lnS,T,r,d,varargin)
%
% Characteristic Function Library of the following models:
%
% Black Scholes
% Heston Stochastic Volatility Model
%

optAlfaCalculation = true;

Error_Message = MException('VerifyInput:InvalidNrOfArguments',...
    'Invalid number of Input arguments');

if strcmp(model,'BlackScholes')
    if nargin == 7
        funobj = @Black_characteristicFn;
    else 
        throw(Error_Message)
    end
elseif strcmp(model,'Heston')
    if nargin == 11
        funobj = @Heston_characteristicFn;
    else 
        throw(Error_Message)
    end
end

fval = feval(funobj,u,lnS,T,r,d,varargin{:});

if optAlfaCalculation == true
    phi = fval;
else
    phi = exp(fval);
end

end


% Explicit Implementation of the characteristic Functions E[exp(iu*lnS_T)]



% Heston
function phi = Heston_characteristicFn(u,lnS,T,r,d,V0,theta,kappa,omega,rho)
    
alfa = -0.5*(u.*u + u*1i);
beta = kappa - rho*omega*u*1i;

omega2 = omega * omega;
gamma = 0.5 * omega2;

D = sqrt(beta .* beta - 4.0 * alfa .* gamma);

bD = beta - D;
eDt = exp(- D * T);


G = bD ./ (beta + D);
B = (bD ./ omega2) .* ((1.0 - eDt) ./ (1.0 - G .* eDt));
psi = (G .* eDt - 1.0) ./(G - 1.0);
A = ((kappa * theta) / (omega2)) * (bD * T - 2.0 * log(psi));


phi = A + B*V0 + 1i*u*(lnS+(r-d)*T);

end
    
