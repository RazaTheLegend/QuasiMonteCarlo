clear; clc;
C = 1;
S = 100; 
K = 1:200; 
T = 1;  
r = 0.0;
d = 0.0;
v = 0.04;  % initial var
theta = 0.04;  
kappa = 0.25;  
sigma = 0.5;  
rho = -0.6; 
some =length(K);
Nsimp = 10000;  %sample paths,

Price1 = ones(1,some);
Price2 = ones(1,some);
Price3 = ones(1,some);
Price4 = ones(1,some);

f = waitbar(0, 'Starting');
for i=1:length(K)
Price1(i) =  RQMC_HestonFullSampling(S_0, r, T, i, V_0, theta, kappa, sigma, rho, Nsimp, CallorPut);
Price2(i) =  QMC_HestonFullSampling(S_0, r, T, i, V_0, theta, kappa, sigma, rho, Nsimp, CallorPut);
Price3(i) =  MyHeston(S_0, r, T, i, V_0, theta, kappa, sigma, rho, Nsimp, 100); %Was chosen as we only used 10000 sample paths
Price4(i) =  FFT_CP('Heston', 18, S_0, i, T, r, d, v, theta, kappa, omega, rho);
    waitbar(i/length(K), f, sprintf('Progress: %d %%', floor(i/length(K)*100)));
    pause(0.1);
end
close(f)

%Plot for option price
figure;
plot(K,Price1(K),'m.',K,Price2(K2),'g.',K,flip(Price3(K)),'b.',K, flip(Price4(K)),'k-')
legend('RQMC','QMC','EMC','FFT')
title("Option prices at different Strikes, Spot 100")
xlabel("K")
ylabel("Option Price")
axis([0 150 -1 100])