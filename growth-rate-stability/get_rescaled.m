%% Rescale to account for growth rate in resilience and recovery
%
% Author: Andrea Mentges, 2020
%
% Input arguments:
%   t = matrix with sampling time points, each column is a different
%   organism
%   N = matrix with biomass from sampling, each column is a different
%   organism
%   r = vector with growth rates of the organisms, each column is a
%   different organism
%   t_dist = time of disturbance
%   varargin = optional input argument, 
%
% Output:
%   t_corr = vector with rescaled time points
%   N_corr = vector with biomass corresponding to rescaled time points

function [t_corr, N_corr] = get_rescaled(t, N, r, t_dist, varargin)

% Rescale by: 
% 1. make time series start at zero (subtract t_dist)
% 2. multiply time by growth rate
% 3. restrict to joint sampling time starting from t_dist

t_corr = (t-t_dist).*repmat(r, size(t,1), 1);
last_joint = min(max(t_corr));
idx = t_corr<=last_joint & t_corr >= 0;
N_corr = N;
N_corr(~idx) = NaN;


end
