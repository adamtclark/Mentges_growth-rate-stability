%% Differential equation for logistic growth function
% where t = time, N = biomass, r = growth rate, and K = carrying capacity

function [Nt] = ode_logistic(t, N, r, K)

% % Constant growth
% Nt = r;

% % Exponential growth
% Nt = r*N;

% % Logistic growth
Nt = r*N*(K-N)/K;

end
