%% Function for solving ODE with a stochastic disturbance regime on top
%
% Time is split into intervals (equidistant or random length).
% At the beginning of each interval, an event happens, which
% influences the starting abundance.
% The integration in each interval takes place at specified, equidistant
% time steps (solver is forced to take these steps).
%
% event_type: The stochastic events are either
% - fixed, each event has the same intensity (I_stoch_mean)
% - random, normally distributed random numbers with mean I_stoch_mean and 
% std I_stoch_std
% - user_defined, the event intensities are user-defined (in I_stoch_mean)
%
% interval_type: the intervals are either
% - fixed, i.e. all the same length, specified by lambda, or
% - random, i.e. of random length, with frequency lambda
% - user_defined, the interval lengths are user-defined (in vector lambda)
%
% Optional input arguments (varagin, variable number of input arguments):
% - 'plot', 'on', will plot default graph
%
% EXAMPLE:
%   solve_stochastic_regime(20, 0.1, 1, 35, 0, 10, 0.03, 'random',...
%       'fixed', 'plot', 'on')
%
% OUTPUT:
% - t, column vector with the time points
% - N, column vector of same size as t containing the respective biomass 
% - event_I, row vector containing for each disturbance event the
% associated disturbance intensity (negative for a reduction of biomass,
% positive for an increase of biomass)
% - event_t, row vector containing the time points at which the
% disturbance events happen

function [t, N, event_I, event_t] = solve_stochastic_regime(tend, r, K,...
    N0, I_stoch_mean, I_stoch_std, lambda, event_type, interval_type,...
    varargin)

% Security checks
assert(length(tend) == 1, 'tend must be one single number')
assert(tend>0, 'tspan must be greater zero')
assert(r>0, 'growth rate r must be greater zero')
assert(K>0, 'capacity K must be greater zero')
assert(N0>0, 'starting abundance N0 must be greater zero')
assert(isnumeric(I_stoch_mean), 'mean disturbance size I_stoch must be numeric')
assert(I_stoch_std>0 || strcmp(event_type, 'user_defined'), 'std of disturbance size I_stoch_std must be greater zero')
assert(all(lambda>0), 'frequency of disturbances must be greater zero')
assert(ischar(event_type), 'event_type must be a string, either fixed or random')
assert(ischar(interval_type), 'interval_type must be a string, either fixed or random')

% By default: no plotting
plotting = 'off';
if any(strcmp(varargin, 'plot'))
    ind  = find(strcmp(varargin, 'plot'));
    plotting = varargin{ind+1};     
end

% No disturbance right at the start of the time series
I_dist = 0; 
t_dist = 0;

% step-width of integration per interval (interval = time between two
% events) i.e. on average, each interval is divided into "step_interval"
% time points
step_interval = 0.1; 

% add t = 0 to vectors
if strcmp(interval_type, 'user_defined')
    lambda = [0 round(lambda*1e8)/1e8 tend];
    I_stoch_mean = [0 I_stoch_mean 0];
end

%% Stochasticity parameters

% Depending on the given event type, assign a function that determines the
% event strength
switch event_type
    case 'fixed'
        f_event = @(t) I_stoch_mean;
    case 'random'
        f_event = @(t) I_stoch_std.*randn(1,1) + I_stoch_mean;
    case 'user_defined'
        f_event = @(t) I_stoch_mean(find(abs(lambda-t)==min(abs(lambda-t))));
end

% Depending on the given interval type, assign a function that determines
% the length of the interval (amount of time passing before next
% disturbance event takes place)
switch interval_type
    case 'fixed'
        f_interval = @(t) lambda;
    case 'random' % draw from exponential distribution with frequency lambda
        %f_interval = @(t) (-1/lambda)*log(rand(1,1)); 
        f_interval = @(t) exprnd(1/lambda, 1, 1); % this and the above 
        % option give the same results!
    case 'user_defined'
        f_interval = @(t) lambda(find(abs(lambda-t)==min(abs(lambda-t)))+1) ...
            - lambda(find(abs(lambda-t)==min(abs(lambda-t))));
end


%% Iterate through intervals

% Starting conditions
tvec = 0;
Nvec = N0;
event_I = [];
event_t = [];

% Keep adding new events until the end of the time series is exceeded
while tvec(end)<=tend
    
    % stochastic event
    t = tvec(end);
    event = f_event(t); 
    N0_interval = Nvec(end) + event; 
    % assert(N0_interval>0, 'starting abundance after event is below zero. reduce event size or standard deviation.')
    
    % length of interval (rounded up to the nearest step size)
    end_interval = ceil(f_interval(t)/step_interval)*step_interval;
    
    % solve differential equation for interval
    % disturbance is assumed to have happened in-between the last time step
    % of last interval and before starting this interval 
    [t_interval, N_interval] = wrap_ode_logistic([0 step_interval end_interval],...
        r, K, N0_interval, I_dist, t_dist, 'plot', 'off');
    if end_interval==step_interval
        t_interval = t_interval([1 end]);
        N_interval = N_interval([1 end]);
    end
    
    % Save event strength and time point for supervision purposes
    event_I = [event_I event];
    event_t = [event_t tvec(end)+step_interval];
    
    % add time and biomass from current interval to output vector
    tvec = [tvec(1:end-1); tvec(end)+t_interval];
    Nvec = [Nvec(1:end-1); N_interval]; 
    
end

% round precision of time to two orders of magnitude lower than step_interval
tvec = round(tvec/(step_interval/100))*(step_interval/100);

% cut the timeseries until desired time point ()
t = tvec(tvec<=tend);
N = Nvec(tvec<=tend);

% Throw error if the time series length is not as expected
assert(length(t)==tend/step_interval+1, 'Weird length of time series.')

%% Plotting
% In case the function is called with the optional input arguments 'plot',
% and 'on', plot the time series

switch plotting
    case {'on'}
        figure('color', 'white', 'position', [-1205,493,426,230])
        p = plot(t, N, '-', 'color', mycolors('blue'));
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
