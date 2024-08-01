%% Get stability measures from timeseries with pulse disturbance
%
% Authors: Andrea Mentges (2020), Adam Clark (2023)
%
% Input arguments:
%   t = vector with sampling time points
%   N = vector with biomass from sampling
%   varargin = optional input argument, either 'paper' (metrics used in the
%   paper), 'simple', or 'Hillebrand2018', which returns measures as based 
%   on the reference below.
%   K = carrying capacity (defaults to 50)
%   
%
% Reference:
% Hillebrand H, Langenheder S, Lebret K, Lindström E, Östman Ö, 
% Striebel M. 2018. Decomposing multiple dimensions of stability in global
% change experiments. Ecology Letters 21:21–30.


function [resilience, recovery, resistance] = get_stability(t, N, K, varargin)
%% Preparations

notnan = ~isnan(t)&~isnan(N);
t = t(notnan);
N = N(notnan);

% find time point of disturbance
idx_dist = find(N==min(N));
t_dist = t(idx_dist);
% t_pre = find(t<t_dist, 1, 'last');
% N_pre = N(t_pre);

N_pre   = N(max(idx_dist-1, 1)); % pre-disturbance biomass
control = N_pre; % N_pre is in equilibrium

% post-disturbance time points
post = t>t_dist;

method = 'paper';
if ~isempty(varargin)
    if strcmp(varargin{:}, 'Hillebrand2018')
        method = 'Hillebrand2018';
    elseif strcmp(varargin{:}, 'paper')
        method = 'paper';
    elseif strcmp(varargin{:}, 'simple')
        method = 'simple';   
    else
        error('undefined method specified.')
    end
end

%% Resilience 

% % My approach: 
% % find the time it takes to re-establish x% of initital biomass
% x = 60;
% idx_t_resil = find(N>=(N_pre*x/100) & t>t_dist, 1, 'first');
% resilience = t(idx_t_resil)-t_dist;
% if isempty(resilience)
%     resilience = NaN;
% end

switch method
    case {'simple', 'Hillebrand2018'}

        % Slope of relative function over time
        x = t(idx_dist:end);
        y = log(N(idx_dist:end)/control);
        % plot(x,y)
        p = polyfit(x,y,1);
        % x1 = linspace(min(x), max(x), 100);
        % % y1 = polyval(p,x1);
        % y1 = p(2) + p(1)*x1;
        % hold on 
        % plot(x1, y1, 'r--')
        resilience = p(1);

        % Note: 
        % maximum resilience: strongly positive
        % minimum resilience: strongly negative, means the function continues to
        % decline
        % at resilience = 0 the function is constant over time
        
    case 'paper'
        
        % get time-to-time slopes
        % slopes = diff(N(post))./diff(t(post));
        
        % resilience = max(slopes);


        % get per capita return rate to equilibrium
        x = N(post) - K;
        xmid = (x(1:(end-1))+x(2:end))/2;
        dnndt = diff(x)./diff(t(post))./xmid;
        
        resilience = max(abs(dnndt));
        
end


%% Recovery

switch method
    case {'simple', 'paper'}
        
        % My approach:
        % percent recovery (relative to pre-disturbance biomass or control)
        recovery = N(end)/N_pre;
    
    case 'Hillebrand2018'

        % LRR compared to control at final sampling point
        recovery = log(N(end)/control);

        % Note: 
        % maximum recovery = 0
        % minimum recovery = minus infinite

    
end


%% Resistance

switch method
    case 'simple'
        
        % My approach:
        % percent persisting after disturbance 
        % (relative to pre-disturbance biomass)
        resistance = 1-(N_pre-N(idx_dist))/N_pre;

        
    case 'Hillebrand2018'
        
        % LRR to pre-disturbance biomass (of control)
        resistance = log(N(idx_dist)/control);

        % Note: 
        % maximum resistence = 0
        % minimum resistence = minus infinite
        
    case 'paper'
        
        % resistance = ratio N/control at first time point after disturbance
        resistance = N(idx_dist)/control;
end

end
