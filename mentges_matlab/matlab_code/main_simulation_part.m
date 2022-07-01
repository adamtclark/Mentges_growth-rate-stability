%% Simulation part: Generating data and plotting figures

clear all
close all

% %%% For publication of code: remove everything followed/preceded by 
% triple percent signs! %%%

% set current folder
cd('/Users/am41xite/Nextcloud/Codes/Matlab_Project_III') %%%

% set figure directory
fdir = '/Users/am41xite/Nextcloud/Figures/2020-09_Stability'; %%%

% set default axes font size
set(groot,'defaultAxesFontSize', 11)
set(groot,'defaultTextFontSize', 12)

% add paths so functions in these folders are found by Matlab
addpath(genpath('/Users/am41xite/Nextcloud/Codes/Matlab_Project_II')) %%%
addpath(genpath('/Users/am41xite/Nextcloud/Codes/Matlab_Functions')) %%%

%% Plotting parameters

% colors
cs = mycolors('lightdarkblue'); % slow-grower
cf = mycolors('yellow'); % fast-grower
cg = mycolors('green2'); % high-capacity ("good") habitat
cb = mycolors('mediumred'); % low-capacity ("bad") habitat
cdist = [.9 .9 .9]; % disturbance bars

% linewidths
ls = 1.5;
lf = 2.5;

% transparency
alpha = 0.5;

%% Global parameters

% growth rates for the two organisms
r_slow = 0.1;
r_fast = 0.5;

% sampling interval for resistance calculation
sf = 3;

% simulation of pulse disturbance (resilience, recovery, resistance)
tend = 20; % sampling duration
tspan = [0 0.01 tend]; % simulation time points
K = 50; % carrying capacity
N0 = K; % starting abundance at the beginning
I_dist = 40; % disturbance intensity
t_dist = 2; % time point of disturbance

% simulation of stochasticity regime (invariability)
tend_stoch = 400; % sampling duration
lambda = 0.03; % frequency of stochastic events
I_stoch_mean = 0; % intensity of stochastic events
I_stoch_std  = 10; % intensity of stochastic events

%% Simulation for a slow- and a fast-grower

%%%% Resilience & Recovery
[t_slow, N_slow] = wrap_ode_logistic(tspan, r_slow, K, N0, I_dist, t_dist, 'plot', 'off');
[t_fast, N_fast] = wrap_ode_logistic(tspan, r_fast, K, N0, I_dist, t_dist, 'plot', 'off');

%%%% Resistance
% sub-sample by sampling frequency
ideally = 1: sf: tspan(end);
[~,idx_resi] = ismember(ideally, t_fast);

%%%% CV
seed = 17; % choose random number generator seed

% % Option 1: simulate a new random time series
% event_type = 'random'; % strength of disturbance events are chosen randomly
% interval_type = 'random'; % length of intervals are chosen randomly
% rng(seed) % set random number seed
% [t_stoch_slow, N_stoch_slow, ~, event_t_slow] = solve_stochastic_regime(tend_stoch, r_slow, K, N0, I_stoch_mean, I_stoch_std, ...
%     lambda, event_type, interval_type, 'plot', 'off');
% rng(seed) % set random number seed
% [t_stoch_fast, N_stoch_fast, event_I_fast, event_t_fast] = solve_stochastic_regime(tend_stoch, r_fast, K, N0, I_stoch_mean, I_stoch_std, ...
%     lambda, event_type, interval_type); 

% Option 2: reproduce the exact time series that is shown in Figure 1 (by
% prescribing the disturbance strenghts and time points)
event_type    = 'user_defined';
interval_type = 'user_defined';
event_I = [-3.95, -15.17,   7.96,    2.86,  -21.68, 17.97, -13.59,  -6.07, 3.70];
event_t = [30,    70.00,  111.00,  125.1,   143.6, 177.9,  231.7, 325.70, 377.0];
rng(seed) % set random number seed
[t_stoch_slow, N_stoch_slow, ~, event_t_slow] = solve_stochastic_regime(tend_stoch, r_slow, K, N0, event_I, NaN, ...
    event_t, event_type, interval_type, 'plot', 'off');
rng(seed) % set random number seed
[t_stoch_fast, N_stoch_fast, ~, event_t_fast] = solve_stochastic_regime(tend_stoch, r_fast, K, N0, event_I, NaN, ...
    event_t, event_type, interval_type);  

%% Figure 1: Effect of growth rate on the four stability metrics

figure('color', 'white', 'position', [155,583,763,165])
subplot(1,3,1)
            ax1 = gca;
            axis(ax1, 'tight')
            ax1.YLim = [0 K*1.1];
            hold on
            dp = plot([t_dist t_dist], ax1.YLim, 'linewidth', 5, 'color', cdist);
            uistack(dp,'bottom')
            hold on
            pf1 = plot(t_fast, N_fast, 'linewidth', lf, 'color', cf);
            ps1 = plot(t_slow, N_slow, 'linewidth', ls, 'color', cs);
            ax1.Layer = 'top';
            ax1.XTick = [];
            ax1.YTick = [];
subplot(1,3,2)
            ax2 = gca;
            axis(ax2, 'tight')
            ax2.YLim = [0 K*1.1];
            hold on
            ps2 = plot(t_slow(idx_resi), N_slow(idx_resi), '--o', 'linewidth', ls, 'color', cs);
            pf  = plot(t_fast, N_fast, '-', 'linewidth', 1, 'color', [.4 .4 .4]);
            pf2 = plot(t_fast(idx_resi), N_fast(idx_resi), '--o', 'linewidth', lf, 'color', cf); 
            dp = plot([t_dist t_dist], ax2.YLim, 'linewidth', 5, 'color', [.9 .9 .9]);
            uistack(dp,'bottom')
            ax2.XTick = [0 10 20];
            ax2.YTick = [];
            ps = plot(t_slow, N_slow, '-', 'linewidth', 1, 'color', [.4 .4 .4]);
            ax2.Layer = 'top';
subplot(1,3,3)
            p = plot(t_stoch_fast, N_stoch_fast, '-', 'color', cf, 'linewidth', lf);
            hold on
            p = plot(t_stoch_slow, N_stoch_slow, '-', 'color', cs, 'linewidth', ls);
            p.Color(4) = alpha;
            ax3 = gca;
            ax3.YLim = [K-I_stoch_std*2.5 K+I_stoch_std*2.5];
            for i = 2:length(event_t_slow)
                dp = plot([event_t_slow(i) event_t_slow(i)], ax3.YLim, 'linewidth', 5, 'color', cdist);
                uistack(dp,'bottom')
            end
            ax3.Layer = 'top';
            ax3.XTick = [0 200 400];
            ax3.YTick = [30 50 70];
            ax3.XLabel.String = 'Time';
            ax3.YLabel.String = 'Biomass';
      
ax1.Position(2) = 0.2;
ax2.Position(2) = 0.2;
ax3.Position(2) = 0.2;            
ax1.Position(3) = 0.17;
ax2.Position(3) = 0.17;
ax3.Position(3) = 0.19;
ax1.Position(4) = 0.62;
ax2.Position(4) = 0.62;
ax3.Position(4) = 0.62;
ax1.Position(1) = 0.1;
ax2.Position(1) = ax1.Position(1)+ax1.Position(3)+0.07;
ax3.Position(1) = ax2.Position(1)+ax2.Position(3)+0.07;

ax1.Box = 'on';
ax2.Box = 'on';
            
l = legend(ax1, [pf1 ps1], 'fast grower','slow grower', 'location', 'eastoutside');
l.Position = [0.61,0.88,0.28,0.06];
l.Box = 'off';

%%%
% export_fig(sprintf('%s/figure1_vary_r_zoom.bmp', fdir), '-r150')

%% Compute stability

method = 'paper';

% Resilience, Recovery
[resil_slow, recov_slow, resis_slow] = get_stability(t_slow, N_slow, method);
[resil_fast, recov_fast, resis_fast] = get_stability(t_fast, N_fast, method);
fprintf('\n\n\t\t Resilience\t Recovery')
fprintf('\nslow\t\t %1.2f\t\t %1.2f\t\t', resil_slow, recov_slow)
fprintf('\nfast\t\t %1.2f\t\t %1.2f', resil_fast, recov_fast)

% Resistance
[~, ~, resis_slow] = get_stability(t_slow, N_slow, method);
[~, ~, resis_fast] = get_stability(t_fast, N_fast, method);
[~, ~, resis_slow2] = get_stability(t_slow(idx_resi), N_slow(idx_resi), method);
[~, ~, resis_fast2] = get_stability(t_fast(idx_resi), N_fast(idx_resi), method);
fprintf('\n\n\t\tResistance')
fprintf('\nslow\t\t %1.2f\t', resis_slow2)
fprintf('\nfast\t\t %1.2f\t', resis_fast2)

% Invariability
CV_slow = std(N_stoch_slow)/mean(N_stoch_slow);
CV_fast = std(N_stoch_fast)/mean(N_stoch_fast);
fprintf('\n\n\t\t1/CV')
fprintf('\nslow\t\t %1.2f', 1/CV_slow)
fprintf('\nfast\t\t %1.2f', 1/CV_fast)

%% Apply partitioning

%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Resilience and recovery
%%%%%%%%%%%%%%%%%%%%%%%%%%%

[t_corr, N_corr] = get_rescaled([t_slow t_fast], [N_slow N_fast], ...
    [r_slow r_fast], t_dist);

t_corr_slow = t_corr(:,1);
t_corr_fast = t_corr(:,2);
N_corr_slow = N_corr(:,1);
N_corr_fast = N_corr(:,2);

[resil_slow, recov_slow] = get_stability(t_corr_slow, N_corr_slow);
[resil_fast, recov_fast] = get_stability(t_corr_fast, N_corr_fast);

fprintf('\n\n\n\t\t Resilience\t Recovery')
fprintf('\nslow\t\t %1.2f\t\t %1.2f\t\t', resil_slow, recov_slow)
fprintf('\nfast\t\t %1.2f\t\t %1.2f\t\t', resil_fast, recov_fast)


%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Resistance (my solution)
%%%%%%%%%%%%%%%%%%%%%%%%%%%

% convert highly-resolved simulation timeseries to "imperfectly sampled" 
% time series (at distinct sampling points)
t_slow_d = t_slow(idx_resi);
N_slow_d = N_slow(idx_resi);
t_fast_d = t_fast(idx_resi);
N_fast_d = N_fast(idx_resi);

% % By the way: the rescaling from above would not help here.
% [t_corr, N_corr] = get_rescaled([t_slow_d t_fast_d], [N_slow_d N_fast_d], ...
%     [r_slow r_fast], t_dist);
% 
% t_corr_slow_d = t_corr(:,1);
% t_corr_fast_d = t_corr(:,1);
% N_corr_slow_d = N_corr(:,1);
% N_corr_fast_d = N_corr(:,1);
% 
% [~,~,resis_slow] = get_stability(t_corr_slow_d, N_corr_slow_d);
% [~,~,resis_slow] = get_stability(t_corr_fast_d, N_corr_fast_d);
% 
% fprintf('\n\n\n\t\t Resistance')
% fprintf('\nslow\t\t %1.2f', resis_slow)
% fprintf('\nfast\t\t %1.2f', resis_fast)

% time passed between first sampling of fast organism and disturbance
first_after_dist = find(t_fast_d>t_dist, 1, 'first');
timestep_fast = t_fast_d(first_after_dist) - t_dist;

% this is the time that should pass between samplings for the slow grower
timestep_slow = timestep_fast*r_fast/r_slow;

% either: use the closest time point to t_dist+t_timestep_slow
% or: interpolate respective time step from slow grower data
x = t_slow_d;
y = N_slow_d;
x2 = t_dist+timestep_slow;
y2 = interp1(x, y, x2, 'linear');

figure()
plot(x, y, 'color', cs)
hold on
plot(t_fast_d, N_fast_d, 'color', cf)
plot(x2, y2, 'bp')
hold off
grid
legend('slow','fast', 'Interpolated Point', 'Location', 'NW')

resis_fast_corr  = 1-(N_fast_d(1)-N_fast_d(first_after_dist))/N_fast_d(1);
resis_slow_corr  = 1-(N_fast_d(1)-y2)/N_fast_d(1);

fprintf('\n\n\t\t Re-scaled Resistance')
fprintf('\nfast\t\t %1.2f', resis_fast_corr)
fprintf('\nslow\t\t %1.2f\n', resis_slow_corr)

%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Invariability
%%%%%%%%%%%%%%%%%%%%%%%%%%%

CV_slow_corr = (std(N_stoch_slow)/mean(N_stoch_slow))*sqrt(r_slow);
CV_fast_corr = (std(N_stoch_fast)/mean(N_stoch_fast))*sqrt(r_fast);
fprintf('\n\t\t 1/CV_corr')
fprintf('\nslow\t\t %1.2f', 1/CV_slow_corr)
fprintf('\nfast\t\t %1.2f', 1/CV_fast_corr)

%% NEW BIT, TO BE TESTED: SUGGESTION BY ADAM, Sep 9th

%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Resistance (NEW solution by Adam)
%%%%%%%%%%%%%%%%%%%%%%%%%%%

% sampling interval for resistance calculation
sf = 2; 

% sub-sample by sampling frequency 
ideally = 2: sf: tspan(end);
[~,idx_resi] = ismember(ideally, t_fast);

% convert highly-resolved simulation timeseries to "imperfectly sampled" 
% time series (at distinct sampling points)
t_slow_d = t_slow(idx_resi);
N_slow_d = N_slow(idx_resi);
t_fast_d = t_fast(idx_resi);
N_fast_d = N_fast(idx_resi);

% time passed between last sampling before disturbance and first sampling 
% after disturbance
t_since_fast = t_fast_d(2) - t_fast_d(1);
t_since_slow = t_slow_d(2) - t_slow_d(1);

% approximate biomass if the disturbance happened right at the last time
% point before disturbance
N_dist_fast = N_fast_d(2)/exp(r_fast*t_since_fast);
N_dist_slow = N_slow_d(2)/exp(r_slow*t_since_slow);

% resulting resistance from that approximated biomass
resis_fast_corr  = 1-(N_fast_d(1)-N_dist_fast)/N_fast_d(1);
resis_slow_corr  = 1-(N_slow_d(1)-N_dist_slow)/N_slow_d(1);

fprintf('\n\n\t\t Re-scaled Resistance (ADAMS METHOD)')
fprintf('\nfast\t\t %1.2f', resis_fast_corr)
fprintf('\nslow\t\t %1.2f\n', resis_slow_corr)

figure('color', 'white', 'position', [155,583,763,165])
            ax1 = gca;
            axis(ax1, 'tight')
            ax1.YLim = [0 K*1.1];
            hold on
            dp = plot([t_dist t_dist], ax1.YLim, 'linewidth', 5, 'color', cdist);
            uistack(dp,'bottom')
            hold on
            pf1 = plot(t_fast, N_fast, 'linewidth', 1, 'color', cf);
            pf1.Color(4) = 0.6;
            ps1 = plot(t_slow, N_slow, 'linewidth', 1, 'color', cs);
            ps1.Color(4) = 0.6;
            pf1d = plot(t_fast_d, N_fast_d,'--o', 'linewidth', 2, 'color', cf);
            ps1d = plot(t_slow_d, N_slow_d,'--o', 'linewidth', 2, 'color', cs);
            pf1d1 = scatter(t_fast_d(2), N_fast_d(2), 35, cf, 'Filled');
            ps1d1 = scatter(t_slow_d(2), N_slow_d(2), 35, cs, 'Filled');
            ax1.Layer = 'top';
            
            sc = scatter(t_dist, N_dist_fast, 35, cf, 'filled', 'marker', '^');
            plot([t_dist t_fast_d(2)], [N_dist_fast N_fast_d(2)], ':', 'color', cf, 'linewidth', 2)
            sc = scatter(t_dist, N_dist_slow, 35, cs, 'filled', 'marker', '^');
            plot([t_dist t_slow_d(2)], [N_dist_slow N_slow_d(2)], ':', 'color', cs, 'linewidth', 2)


%% Figure 3: Illustrate resilience & recovery rescaling


[t1, N1] = wrap_ode_logistic(tspan, 1, K, N0, I_dist, t_dist, 'plot', 'off');

figure('color', 'white', 'position', [155,583,763,165])
subplot(1,3,1)
            ax1 = gca;
            axis(ax1, 'tight')
            ax1.YLim = [0 K*1.1];
            hold on
            dp = plot([t_dist t_dist], ax1.YLim, 'linewidth', 5, 'color', cdist);
            uistack(dp,'bottom')
            hold on
            pf1 = plot(t_fast, N_fast, 'linewidth', lf+3, 'color', cf);
            ps1 = plot(t_slow, N_slow, 'linewidth', ls, 'color', cs);
            ax1.Layer = 'top';
subplot(1,3,2)
            ax2 = gca;
            axis(ax2, 'tight')
            ax2.YLim = [0 K*1.1];
            hold on
            pa = patch([0 1.75 1.75 0],[ax2.YLim(2) ax2.YLim(2) 0 0], ...
                mycolors('paleblue'), 'FaceAlpha', 0.5, 'LineStyle', 'none');
            hh = hatchfill(pa, 'single', 10, 5, 'none');
            hh.Color = [.6 .6 .6];
            uistack(pa,'bottom')
            dp = plot([0 0], ax1.YLim, 'linewidth', 5, 'color', [.9 .9 .9]);
            uistack(dp,'bottom')
            hold on
            pf2 = plot((t_fast-t_dist).*repmat(r_fast, size(t_fast,1), 1), N_fast, 'linewidth', lf+3, 'color', cf);
            ps2 = plot((t_slow-t_dist).*repmat(r_slow, size(t_slow,1), 1), N_slow, 'linewidth', ls, 'color', cs);
            ax1.Layer = 'top';   
%             plot(t1-t_dist, N1, 'linewidth', lf+3, 'color', 'red', 'linestyle', ':');
            
ax1.XLabel.String = 'Time';
ax2.XLabel.String = 'Rescaled time';
ax1.YLabel.String = 'Biomass';
      
ax1.Position(2) = 0.24;
ax2.Position(2) = 0.24;           
ax1.Position(3) = 0.17;
ax2.Position(3) = 0.17;
ax1.Position(4) = 0.62;
ax2.Position(4) = 0.62;
ax1.Position(1) = 0.1;
ax2.Position(1) = ax1.Position(1)+ax1.Position(3)+0.07;

ax1.Box = 'on';
ax2.Box = 'on';

annotation('textbox', [0.095, 0.99, 0, 0], 'string', 'a', 'FontWeight', 'bold', 'Fontsize', 12)
annotation('textbox', [0.335, 0.99, 0, 0], 'string', 'b', 'FontWeight', 'bold', 'Fontsize', 12)

annotation('textbox', [0.095+0.02, 0.978, 0.6, 0], 'string', 'Realized', 'FontWeight', 'normal', 'Fontsize', 10.5, 'LineStyle', 'none')
annotation('textbox', [0.335+0.02, 0.978, 0.6, 0], 'string', 'Intrinsic', 'FontWeight', 'normal', 'Fontsize', 10.5, 'LineStyle', 'none')
        
l = legend(ax1, [pf1 ps1], 'fast-grower','slow-grower', 'subset to use', 'location', 'eastoutside');
l.Position = [0.455347313237222,0.708484848484848,0.280000000000000,0.184848484848485];
l.Box = 'off';

%%%
% export_fig(sprintf('%s/illustrate_resilience_rescaling.bmp', fdir), '-r150')

%% Figure 4: Illustrate resistance rescaling

idx_fast_new = ismember(t_slow, [1 3.5:1.5:19.5]);
idx_slow_new = ismember(t_slow, [1 10.25 18]);

t_fast_new = t_fast(idx_fast_new);
N_fast_new = N_fast(idx_fast_new);
t_slow_new = t_slow(idx_slow_new);
N_slow_new = N_slow(idx_slow_new);
t_fslow_new = t_slow(idx_fast_new);
N_fslow_new = N_slow(idx_fast_new);

figure('color', 'white', 'position', [155,583,763,165])
subplot(1,3,1)
            ax1 = gca;
            axis(ax1, 'tight')
            ax1.YLim = [0 K*1.1];
            hold on
            dp = plot([t_dist t_dist], ax1.YLim, 'linewidth', 5, 'color', cdist);
            uistack(dp,'bottom')
            hold on
            pf1 = plot(t_fast, N_fast, 'linewidth', 1, 'color', cf);
            pf1.Color(4) = 0.6;
            ps1 = plot(t_slow, N_slow, 'linewidth', 1, 'color', cs);
            ps1.Color(4) = 0.6;
            pf1d = plot(t_fast_new, N_fast_new,'--o', 'linewidth', 2, 'color', cf);
            ps1d = plot(t_fslow_new, N_fslow_new,'--o', 'linewidth', 2, 'color', cs);
            pf1d1 = scatter(t_fast_new(2), N_fast_new(2), 35, cf, 'Filled');
            ps1d1 = scatter(t_fslow_new(2), N_fslow_new(2), 35, cs, 'Filled');
            ax1.Layer = 'top';
subplot(1,3,2)
            ax2 = gca;
            axis(ax2, 'tight')
            ax2.YLim = [0 K*1.1];
            hold on
            dp = plot([t_dist t_dist], ax1.YLim, 'linewidth', 5, 'color', cdist);
            uistack(dp,'bottom')
            hold on
            pf2 = plot(t_fast, N_fast, 'linewidth', 1, 'color', cf);
            pf2.Color(4) = 0.6;
            ps2 = plot(t_slow, N_slow, 'linewidth', 1, 'color', cs);
            ps2.Color(4) = 0.6;
            pf2d = plot(t_fast_new, N_fast_new,'--o', 'linewidth', 2, 'color', cf);
            ps2d = plot(t_slow_new, N_slow_new,'--o', 'linewidth', 2, 'color', cs);
            pf2d1 = scatter(t_fast_new(2), N_fast_new(2), 35, cf, 'Filled');
            ps2d1 = scatter(t_slow_new(2), N_slow_new(2), 35, cs, 'Filled');
            
            ax2.Layer = 'top';            
            
ax1.XLabel.String = 'Time';
ax2.XLabel.String = 'Time';
ax1.YLabel.String = 'Biomass';
      
ax1.Position(2) = 0.24;
ax2.Position(2) = 0.24;           
ax1.Position(3) = 0.17;
ax2.Position(3) = 0.17;
ax1.Position(4) = 0.62;
ax2.Position(4) = 0.62;
ax1.Position(1) = 0.1;
ax2.Position(1) = ax1.Position(1)+ax1.Position(3)+0.09;

ax1.Box = 'on';
ax2.Box = 'on';

% vertical lines indicating resistance in subplot a
xf = [0.28 0.28];
yf = [0.435 0.805];
xs = [0.29 0.29];
ys = [0.34 0.805];
annotation('line', xf, yf, 'Color', cf, 'Linestyle', '-')
annotation('line', [xf(1)-0.003 xf(1)+0.003], [yf(1) yf(1)], 'Color', cf, 'Linestyle', '-')
annotation('line', [xf(1)-0.003 xf(1)+0.003], [yf(2) yf(2)], 'Color', cf, 'Linestyle', '-')
annotation('line', xs, ys, 'Color', cs, 'Linestyle', '-')
annotation('line', [xs(1)-0.003 xs(1)+0.003], [ys(1) ys(1)], 'Color', cs, 'Linestyle', '-')
annotation('line', [xs(1)-0.003 xs(1)+0.003], [ys(2) ys(2)], 'Color', cs, 'Linestyle', '-')

% vertical lines indicating resistance in subplot b
xf = [0.54 0.54];
yf = [0.435 0.805];
xs = [0.55 0.55];
ys = [0.435 0.805];
annotation('line', xf, yf, 'Color', cf, 'Linestyle', '-')
annotation('line', [xf(1)-0.003 xf(1)+0.003], [yf(1) yf(1)], 'Color', cf, 'Linestyle', '-')
annotation('line', [xf(1)-0.003 xf(1)+0.003], [yf(2) yf(2)], 'Color', cf, 'Linestyle', '-')
annotation('line', xs, ys, 'Color', cs, 'Linestyle', '-')
annotation('line', [xs(1)-0.003 xs(1)+0.003], [ys(1) ys(1)], 'Color', cs, 'Linestyle', '-')
annotation('line', [xs(1)-0.003 xs(1)+0.003], [ys(2) ys(2)], 'Color', cs, 'Linestyle', '-')

% subplot labels (a, b)
annotation('textbox', [0.095, 0.99, 0, 0], 'string', 'a', 'FontWeight', 'bold', 'Fontsize', 12)
annotation('textbox', [0.355, 0.99, 0, 0], 'string', 'b', 'FontWeight', 'bold', 'Fontsize', 12)

% subplot titles
annotation('textbox', [0.095+0.02, 0.975, 0.6, 0], 'string', 'Realized resistance', 'FontWeight', 'normal', 'Fontsize', 10.5, 'LineStyle', 'none')
annotation('textbox', [0.355+0.02, 0.975, 0.6, 0], 'string', 'Intrinsic resistance', 'FontWeight', 'normal', 'Fontsize', 10.5, 'LineStyle', 'none')
            
l = legend(ax1, [pf1 ps1], 'fast-grower','slow-grower', 'subset to use', 'location', 'eastoutside');
l.Position = [0.5,0.708484848484848,0.280000000000000,0.184848484848485];
l.Box = 'off';

%%%
% export_fig(sprintf('%s/illustrate_resistance_rescaling.bmp', fdir), '-r150')

%% Figure 2: Illustrate CV correction

% % Option A) Generate new simulation data with following code:
%
% nvar = 1000;
% growth_rates = linspace(0.05, 2, nvar);
% 
% vars = NaN(1,nvar);
% means = NaN(1,nvar);
% CVs = NaN(1,nvar);
% 
% rng(1)
% for i = 1:nvar
%     i
%     r = growth_rates(i);
%     [t_stoch, N_stoch, event_I, event_t] = solve_stochastic_regime(tend_stoch, r, K, N0, I_stoch_mean, I_stoch_std, ...
%         lambda, event_type, interval_type, 'plot', 'off');
%     
%     vars(i)  = var(N_stoch);
%     means(i) = mean(N_stoch);
%     CVs(i)   = std(N_stoch)/mean(N_stoch);
% end
% 
% % save('workspace_CV_correction_2020-04-21_nvars1000')

% Option B: Load simulation data
load('workspace_CV_correction_2020-04-21')
% load('workspace_CV_correction_2020-04-21_tend-double') %%%
% load('workspace_CV_correction_2020-04-21_lambda-double') %%%
% load('workspace_CV_correction_2020-04-21_nvars1000') %%%

invariability = 1./CVs;

figure('color', 'white', 'position', [339,667,953,146])
subplot(1,4,1)
p = plot(growth_rates, invariability, '.-', 'color', [.8 .8 .8],...
    'MarkerSize', 8, 'MarkerEdgeColor', [.5 .5 .5]);
ax1 = gca;
hold on
l = plot(growth_rates, 1./sqrt(((I_stoch_mean^2)+(I_stoch_std^2))*...
    (lambda./(2*growth_rates))./(means.^2)),...
    '-', 'linewidth', 3, 'Color', mycolors('redorange'));
l.Color(4) = 0.7;
xlabel('Growth rate'), ylabel('Temp. stab._{real.}')
axis tight

CVs_corr = CVs.*sqrt(growth_rates);
invariability_corr = (1./CVs).*(1./sqrt(growth_rates));

subplot(1,4,2)
p = plot(growth_rates, invariability_corr, '.-', 'color', [.8 .8 .8],...
    'MarkerSize', 8, 'MarkerEdgeColor', [.5 .5 .5]);
ax2 = gca;
hold on
l = plot(growth_rates, 1./sqrt(((I_stoch_mean^2)+(I_stoch_std^2)) *...
    lambda./(2*(means.^2))),...
    '-', 'linewidth', 3, 'Color', mycolors('redorange'));
l.Color(4) = 0.7;
xlabel('Growth rate'), ylabel('Temp. stab._{intr.}')
axis tight

linkaxes([ax1 ax2], 'y')
ax2.YLim = [0 100];

ax1.Position(3) = 0.14;
ax2.Position(3) = 0.14;
ax1.Position(4) = 0.65;
ax2.Position(4) = 0.65;

% Subpanel labels (a, b)
annotation('textbox', [0.127, 0.999, 0, 0], 'string', 'a', ...
    'FontWeight', 'bold', 'Fontsize', 11)
annotation('textbox', [0.333, 0.999, 0, 0], 'string', 'b', ...
    'FontWeight', 'bold', 'Fontsize', 11)

% Subpanel titles
annotation('textbox', [0.145, 0.999, 0.6, 0], 'string',...
    'Realized', 'FontWeight', 'normal', 'Fontsize', 11,...
    'LineStyle', 'none')
annotation('textbox', [0.35, 0.999, 0.6, 0], 'string',...
    'Intrinsic', 'FontWeight', 'normal', 'Fontsize', 11,...
    'LineStyle', 'none')

l = legend('Stochastic simulations', 'Analytical solution',...
    'Box', 'off');
l.Position = [0.499475341028332,0.674657534246575,0.1594998,0.20890441];

std(CVs)

%%%
% export_fig(sprintf('%s/illustrate_CV_solution.bmp', fdir), '-r150')


%% Supplement: Conceptual simulation of effect of K

%%%%%%%%%
% Resilience & Recovery
%%%%%%%%%
r = r_slow;
K_good = 80;
K_bad = 50;
N0_good = K_good;
N0_bad = K_bad;

[t_good, N_good] = wrap_ode_logistic(tspan, r, K_good, N0_good,...
    I_dist, t_dist, 'plot', 'off');
[t_bad,  N_bad]  = wrap_ode_logistic(tspan, r, K_bad,  N0_bad,...
    I_dist, t_dist, 'plot', 'off');

%%%%%%%%%
% Resistance
%%%%%%%%%
% sub-sample by sampling frequency
ideally = 1: sf: tspan(end);
[~,idx_resi] = ismember(ideally, t_fast);

%%%%%%%%%
% CV
%%%%%%%%%
rng(seed) % set random number seed
[t_stoch_good, N_stoch_good, ~, event_t_good] = solve_stochastic_regime(...
    tend_stoch, r, K_good, N0_good, I_stoch_mean, I_stoch_std, ...
    lambda, event_type, interval_type);
rng(seed) % set random number seed
[t_stoch_bad , N_stoch_bad , ~, event_t_bad ] = solve_stochastic_regime(...
    tend_stoch, r, K_bad , N0_bad , I_stoch_mean, I_stoch_std, ...
    lambda, event_type, interval_type);
      
%% Figure S3: Plot effect of K on stability

figure('color', 'white', 'position', [155,583,763,165])
subplot(1,3,1)
            ax1 = gca;
            axis(ax1, 'tight')
            ax1.YLim = [0 K_good*1.1];
            ax1.YLabel.String = 'Biomass';
            hold on
            dp = plot([t_dist t_dist], ax1.YLim, 'linewidth', 5, 'color', [.9 .9 .9]);
            uistack(dp,'bottom')
            text(t_dist+3, ax1.YLim(2)*1.14, 'Disturbance', 'HorizontalAlignment', 'Center', ...
               'Fontsize', 10, 'color', [.7 .7 .7]) 
            text(t_dist, ax1.YLim(2)*1.04, '▼', 'HorizontalAlignment', 'Center', ...
               'Fontsize', 8, 'color', [.7 .7 .7]) 
            hold on
            pf1 = plot(t_good, N_good, 'linewidth', lf, 'color', cg);
            ps1 = plot(t_bad, N_bad, 'linewidth', ls, 'color', cb);
            ax1.Layer = 'top';
            ax1.XTick = [0 10 20];
subplot(1,3,2)
            ax2 = gca;
            axis(ax2, 'tight')
            ax2.YLim = [0 K_good*1.1];
            hold on
            ps2 = plot(t_bad(idx_resi), N_bad(idx_resi), '--o', 'linewidth', ls, 'color', cb);
            pf  = plot(t_good, N_good, '-', 'linewidth', 1, 'color', [.4 .4 .4]);
            pf2 = plot(t_good(idx_resi), N_good(idx_resi), '--o', 'linewidth', lf, 'color', cg); 
            dp = plot([t_dist t_dist], ax2.YLim, 'linewidth', 5, 'color', [.9 .9 .9]);
            uistack(dp,'bottom')
            ax2.XLabel.String = 'Time';
            ax2.XTick = [0 10 20];
            ax2.YTick = [];
            ps = plot(t_bad, N_bad, '-', 'linewidth', 1, 'color', [.4 .4 .4]);
            ax2.Layer = 'top';
subplot(1,3,3)
            p = plot(t_stoch_good, N_stoch_good, '-', 'color', cg, 'linewidth', lf);
            hold on
            p = plot(t_stoch_bad, N_stoch_bad, '-', 'color', cb, 'linewidth', ls);
            ax3 = gca;
            ax3.YLim = [K_bad-I_stoch_std*3 K_good+I_stoch_std*2.6];
            for i=1:length(event_t_bad)
                dp = plot([event_t_bad(i) event_t_bad(i)], ax3.YLim, 'linewidth', 5, 'color', [.93 .93 .93]);
                uistack(dp,'bottom')
            end
            ax3.Layer = 'top';
            ax3.XTick = [0 200 400];
            ax3.YTick = [40 80];
      
ax1.Position(2) = 0.202;
ax2.Position(2) = 0.202;
ax3.Position(2) = 0.202;            
ax1.Position(3) = 0.19;
ax2.Position(3) = 0.19;
ax3.Position(3) = 0.19;
ax1.Position(4) = 0.62;
ax2.Position(4) = 0.62;
ax3.Position(4) = 0.62;
ax1.Position(1) = 0.1;
ax2.Position(1) = ax1.Position(1)+ax1.Position(3)+0.07;
ax3.Position(1) = ax2.Position(1)+ax2.Position(3)+0.07;

ax1.Box = 'on';
ax2.Box = 'on';

% subpanels (abc)
annotation('textbox', [0.095, 0.93, 0, 0], 'string', 'a', 'FontWeight', 'bold', 'Fontsize', 11)
annotation('textbox', [0.355, 0.93, 0, 0], 'string', 'b', 'FontWeight', 'bold', 'Fontsize', 11)
annotation('textbox', [0.615, 0.93, 0, 0], 'string', 'c', 'FontWeight', 'bold', 'Fontsize', 11)

% annotations (high resilience, low resilience, ..)
annotation('textbox', [0.1, 0.72, 0.1, 0], 'string', sprintf('high\nresilience'),...
    'Fontsize', 9, 'Color', cg, 'LineStyle', 'none','HorizontalAlignment', 'right');
annotation('textbox', [0.19, 0.83, 0.1, 0], 'string', sprintf('high\nrecovery'),...
    'Fontsize', 9, 'Color', cg, 'LineStyle', 'none','HorizontalAlignment', 'right');
annotation('textbox', [0.025, 0.43, 0.2, 0], 'string', sprintf('low resilience'),...
    'Fontsize', 9, 'Color', cb, 'LineStyle', 'none', 'HorizontalAlignment', 'right');
annotation('textbox', [0.085, 0.38, 0.2, 0], 'string', sprintf('low\nrecovery'),...
    'Fontsize', 9, 'Color', cb, 'LineStyle', 'none', 'HorizontalAlignment', 'right');
annotation('textbox', [0.42, 0.5, 0.35, 0], 'string', 'high resistance',...
    'Fontsize', 9, 'Color', cg, 'LineStyle', 'none');
annotation('textbox', [0.42, 0.295, 0.35, 0], 'string', 'low resistance',...
    'Fontsize', 9, 'Color', cb, 'LineStyle', 'none');
annotation('textbox', [0.67, 0.84, 0.35, 0], 'string', 'high temporal stability',...
    'Fontsize', 9, 'Color', cg, 'LineStyle', 'none');
annotation('textbox', [0.62, 0.3, 0.35, 0], 'string', 'low temporal stability',...
    'Fontsize', 9, 'Color', cb, 'LineStyle', 'none');
            
l = legend(ax1, [pf1 ps1], 'high capacity','low capacity', 'location', 'eastoutside');
l.Position = [0.61,0.88,0.28,0.06];
l.Box = 'off';

%%%
% export_fig(sprintf('%s/figure2_vary_K.bmp', fdir), '-r150')

%% Table S1: Compute stability for simulation of K

method = 'paper';

% Resilience, Recovery
[resil_bad, recov_bad, resis_bad] = get_stability(t_bad, N_bad, method);
[resil_good, recov_good, resis_good] = get_stability(t_good, N_good, method);
fprintf('\n\n\t\t Resilience\t Recovery')
fprintf('\nbad\t\t %1.2f\t\t %1.2f\t\t', resil_bad, recov_bad)
fprintf('\ngood\t\t %1.2f\t\t %1.2f', resil_good, recov_good)

% Resistance
[~, ~, resis_bad] = get_stability(t_bad, N_bad, method);
[~, ~, resis_good] = get_stability(t_good, N_good, method);
[~, ~, resis_bad2] = get_stability(t_bad(idx_resi), N_bad(idx_resi), method);
[~, ~, resis_good2] = get_stability(t_good(idx_resi), N_good(idx_resi), method);
fprintf('\n\n\t\tResistance')
fprintf('\nbad\t\t %1.2f\t', resis_bad2)
fprintf('\ngood\t\t %1.2f\t', resis_good2)

% Invariability
CV_bad = std(N_stoch_bad)/mean(N_stoch_bad);
CV_good = std(N_stoch_good)/mean(N_stoch_good);
fprintf('\n\n\t\t1/CV')
fprintf('\nbad\t\t %1.2f', 1/CV_bad)
fprintf('\ngood\t\t %1.2f', 1/CV_good)



