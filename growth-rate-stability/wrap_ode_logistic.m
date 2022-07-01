%% Wrap function for solving ODE of logistic growth model
% Assigns input arguments and calls the ODE solver.
%
% INPUT:
% tspan, vector containing either [start, end] of time series or 
%   [start, stepwidth, end]. In the first case, if no stepwidth is 
%   provided, the solver chooses a (non-equidistant) stepwidth by itself.
%   In the second case, the solver is forced to evaluate the equation at
%   the specified (equidistant) time steps.
% r, growth rate in logistic growth model
% K, carrying capacity in logistic growth model
% N0, starting biomass at the beginning of the simulation
% I_dist, disturbance intensity, i.e. the amount of biomass by which N is
%   reduced at the time of disturbance
% t_dist, disturbance time point
%
% Optional input arguments (varargin): 
%   'plot', 'on': a figure showing the time series is plotted
%
% OUTPUT: 
% t, column vector with time points
% N, column vector with associated biomass
%
% EXAMPLES: 
% wrap_ode_logistic([0 20], 0.1, 50, 45, 40, 2, 'plot', 'on')

function [t, N] = wrap_ode_logistic(tspan, r, K, N0, I_dist, t_dist, varargin)
%% Check and assign input

% assert(t_dist>0, 't_dist must be positive. Set I_dist=0 if you dont want disturbance.')
assert(K>0, 'K must be positive.')
assert(tspan(1) == 0, 'Results are weird if t>0!')


%% Call ODE solver 

% ODE solver options
options = odeset('NonNegative', 1:size(N0,1));

% Call the solver
if length(tspan)==2 % let solver choose timesteps (non-interpolated)
    if t_dist > 0
        assert(t_dist<tspan(2), 'The disturbance must occur before the end of the simulation.')
        
        % solve until disturbance ("pre")
        [t_pre, N_pre]   = ode45(@ode_logistic, [tspan(1) t_dist], N0, options, r, K);
        
        % solve from disturbance to end ("post")
        [t_post, N_post] = ode45(@ode_logistic, [t_dist+0.01 tspan(2)], N_pre(end)-I_dist, options, r, K); 
        
        % concatenate
        t = [t_pre; t_post];
        N = [N_pre; N_post];
        
    else
        % no disturbance
        [t, N]   = ode45(@ode_logistic, tspan, N0, options, r, K);
    end

elseif length(tspan)==3 % Force the solver to evaluate at specified time steps
    
    if t_dist > 0 
        % solve until disturbance ("pre")
        [t_pre, N_pre]   = ode45_(@ode_logistic, tspan(1):tspan(2):t_dist, N0, options, r, K);

        % solve from disturbance to end ("post")
        [t_post, N_post] = ode45_(@ode_logistic, t_dist+0.01:tspan(2):tspan(3), N_pre(end)-I_dist, options, r, K); 

        % concatenate
        t = [t_pre; t_post];
        N = [N_pre; N_post];
        
    else
        [t, N]   = ode45_(@ode_logistic, tspan(1):tspan(2):tspan(3), N0-I_dist, options, r, K);
        
    end
    
else
    error('Specify tspan as either [tstart tend] or [tstart step tend]!')
end


%% Plotting

if any(strcmp(varargin, 'plot'))
    ind  = find(strcmp(varargin, 'plot'));
    switch varargin{ind+1}
        case {'on'}
            figure('color', 'white', 'position', [-1205,493,426,230])
            plot(t, N, '.-', 'linewidth', 1);
            ax = gca;
            axis(ax, 'tight'), xlabel('time')
            ax.YLim = [0 K*1.1];
            ax.YLabel.String = 'N';
            hold on
            if I_dist>0
                dp = plot([t_dist t_dist], ax.YLim, 'color', [.7 .7 .7]);
                uistack(dp,'bottom')
            end
            set(findall(gcf,'-property','FontSize'),'FontSize',13)
    end
end


end

