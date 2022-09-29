%
% Simple function to calculate the cashflow seen in RQMC
%

function optval = Call_Put(S,K,C)
% S = Nsimp x Nt matrix of simulated prices
% K = Strike price
% C = 1 -> Call; C = 0 -> Put

    Nt = size(S,2);         % Number of Discretization Steps    
    if(C==1)
        optval = mean(max(S(:,Nt)-K,0));
    else
        optval = mean(max(K-S(:,Nt),0));
    end
end