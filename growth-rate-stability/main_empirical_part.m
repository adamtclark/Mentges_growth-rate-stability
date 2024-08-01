%% Empirical part: loading data and plotting figures
%
% Authors: Andrea Mentges (2019) and Adam Clark (2023)

clear all

fdir = mydirectory; % paste user directory here;
addpath(genpath('MYDIRECTORY')) % add folder with code to path 

%% Load data

% load Hillebrand & Kunze (2020) meta-data 
[M, Table_S1] = load_hillebrand2020();

% adds growth rate estimates for each organism group to M
[M, Table_S3] = load_growth_rates(M);

% adds stability estimates (based on the log-response ratios available at
% datadryad.org/stash/dataset/doi:10.5061/dryad.cz8w9gj09)
[M, Table_S2] = load_stability_data(M);


%% Figure 5: Boxplots observed stability

figure('color', 'white', 'position', [45,555,1146,380])
subplot(1,4,1)
[ax1] = fplot_boxplot_log(M.ratio_recov(~isnan(M.growth)), M.growth(~isnan(M.growth)));
ax1.YLabel.String = 'Recovery_{realiz.}';
ax1.YLim = [-0.5 3.5];
lines =  findobj(gca,'Type','line');
l = legend(ax1, [lines(1) lines(2)], 'Median', 'Mean');
l.Position = [0.05 0.77 0.07 0.09];
l.Box = 'off';

subplot(1,4,2)
[ax2] = fplot_boxplot_log(M.ratio_resil(~isnan(M.growth)), M.growth(~isnan(M.growth)));
ax2.YLabel.String = 'Resilience_{realiz.}';
ax2.YLim = [-0.2 1.3];

subplot(1,4,3)
[ax3] = fplot_boxplot_log(M.ratio_resist(~isnan(M.growth)), M.growth(~isnan(M.growth)));
ax3.YLabel.String = 'Resistance_{realiz.}';
ax3.YLim = [-0.1 2.5];

subplot(1,4,4)
[ax4] = fplot_boxplot_log(M.ratio_temp_stab(~isnan(M.growth)), M.growth(~isnan(M.growth)));
ax4.YLabel.String = 'Temporal stability_{realiz.}';
ax4.YLim = [0.3 50];
ax4.YScale = 'log';

annotation('textbox', [0.045, 0.77, 0, 0], 'string', 'a', 'FontWeight', 'bold')
annotation('textbox', [0.265, 0.77, 0, 0], 'string', 'b', 'FontWeight', 'bold')
annotation('textbox', [0.485, 0.77, 0, 0], 'string', 'c', 'FontWeight', 'bold')
annotation('textbox', [0.705, 0.77, 0, 0], 'string', 'd', 'FontWeight', 'bold')

set(findall(gcf,'-property','FontSize'),'FontSize',16)

n_recov = sum(~isnan(M.ratio_recov(~isnan(M.growth))));
n_resil = sum(~isnan(M.ratio_resil(~isnan(M.growth))));
n_resist = sum(~isnan(M.ratio_resist(~isnan(M.growth))));
n_temp_stab = sum(~isnan(M.ratio_temp_stab(~isnan(M.growth))));

annotation('textbox', [0.155, 0.5, 0.2, 0.2], 'string', ...
    sprintf('n = %d', n_recov), 'Fontsize', 12, 'linestyle', 'none')
annotation('textbox', [0.375, 0.5, 0.2, 0.2], 'string', ...
    sprintf('n = %d', n_resil), 'Fontsize', 12, 'linestyle', 'none')
annotation('textbox', [0.375+0.22, 0.5, 0.2, 0.2], 'string', ...
    sprintf('n = %d', n_resist), 'Fontsize', 12, 'linestyle', 'none')
annotation('textbox', [0.375+0.22*2, 0.5, 0.2, 0.2], 'string', ...
    sprintf('n = %d', n_temp_stab), 'Fontsize', 12, 'linestyle', 'none')

% add organism group labels
xpos = unique(M.growth(~isnan(M.growth)));
for o = 1:length(xpos)
    grpname = unique(M.organism(M.growth==xpos(o)));
    ht = text(ax4, xpos(o), ax4.YLim(2)*1.2, grpname, 'FontSize', 12);
    set(ht,'Rotation',45)
end

ax1.Position(1) = 0.05;
ax2.Position(1) = ax1.Position(1) + 0.22;
ax3.Position(1) = ax1.Position(1) + 0.22*2;
ax4.Position(1) = ax1.Position(1) + 0.22*3;
ax1.Position(2) = 0.25;
ax2.Position(2) = 0.25;
ax3.Position(2) = 0.25;
ax4.Position(2) = 0.25;
ax1.Position(3) = 0.15;
ax2.Position(3) = 0.15;
ax3.Position(3) = 0.15;
ax4.Position(3) = 0.15;
ax1.Position(4) = 0.45;
ax2.Position(4) = 0.45;
ax3.Position(4) = 0.45;
ax4.Position(4) = 0.45;

ax1.XLabel.Position(2) = -30;
ax2.XLabel.Position(2) = -30;
ax3.XLabel.Position(2) = -30;
ax4.XLabel.Position(2) = -30;


%% Figure S1: Boxplots - Four different resilience metrics

figure('color', 'white', 'position', [45,555,1146,380])
subplot(1,4,1)
[ax1] = fplot_boxplot_log(M.resil_rma_day(~isnan(M.growth)), M.growth(~isnan(M.growth)));
ax1.YLabel.String = [{'\bfResilience'}; {'\rmas in Hillebrand & Kunze'}];
ax1.YLim = [-1 2];
lines =  findobj(gca,'Type','line');
l = legend(ax1, [lines(1) lines(2)], 'Median', 'Mean');
l.Position = [0.05 0.77 0.07 0.09];
l.Box = 'off';

subplot(1,4,2)
[ax2] = fplot_boxplot_log(M.ratio_resil_regression_log(~isnan(M.growth)), M.growth(~isnan(M.growth)));
ax2.YLabel.String = [{'\rmas regression over log(time)'}];
ax2.YLim = [-2 1];

subplot(1,4,3)
[ax3] = fplot_boxplot_log(M.ratio_resil_regression_lin(~isnan(M.growth)), M.growth(~isnan(M.growth)));
ax3.YLabel.String = [{'\rmas regression over time'}];
ax3.YLim = [-0.03 0.03];

subplot(1,4,4)
[ax4] = fplot_boxplot_log(M.ratio_resil(~isnan(M.growth)), M.growth(~isnan(M.growth)));
ax4.YLabel.String = [{'\rmas maximum slope'}];
ax4.YLim = [-0.05 0.3];

annotation('textbox', [0.045, 0.77, 0, 0], 'string', 'a', 'FontWeight', 'bold')
annotation('textbox', [0.265, 0.77, 0, 0], 'string', 'b', 'FontWeight', 'bold')
annotation('textbox', [0.485, 0.77, 0, 0], 'string', 'c', 'FontWeight', 'bold')
annotation('textbox', [0.705, 0.77, 0, 0], 'string', 'd', 'FontWeight', 'bold')

set(findall(gcf,'-property','FontSize'),'FontSize',16)

% add organism group labels
xpos = unique(M.growth(~isnan(M.growth)));
for o = 1:length(xpos)
    grpname = unique(M.organism(M.growth==xpos(o)));
    ht = text(ax4, xpos(o), ax4.YLim(2)*1.1, grpname, 'FontSize', 12);
    set(ht,'Rotation',45)
end

ax1.Position(1) = 0.06;
ax2.Position(1) = ax1.Position(1) + 0.22;
ax3.Position(1) = ax1.Position(1) + 0.22*2;
ax4.Position(1) = ax1.Position(1) + 0.22*3;
ax1.Position(2) = 0.25;
ax2.Position(2) = 0.25;
ax3.Position(2) = 0.25;
ax4.Position(2) = 0.25;
ax1.Position(3) = 0.15;
ax2.Position(3) = 0.15;
ax3.Position(3) = 0.15;
ax4.Position(3) = 0.15;
ax1.Position(4) = 0.45;
ax2.Position(4) = 0.45;
ax3.Position(4) = 0.45;
ax4.Position(4) = 0.45;

ax1.XLabel.Position(2) = -30;
ax2.XLabel.Position(2) = -30;
ax3.XLabel.Position(2) = -30;
ax4.XLabel.Position(2) = -30;


%% Figure 6: Boxplots - Corrected stability

figure('color', 'white', 'position', [45,555,1146,380])
subplot(1,4,1)
[ax1] = fplot_boxplot_log(M.ratio_recov_corr(~isnan(M.growth)), M.growth(~isnan(M.growth)));
ax1.YLabel.String = 'Recovery_{intr.}';
ax1.YLim = [-0.5 3.5];
lines =  findobj(gca,'Type','line');
l = legend(ax1, [lines(1) lines(2)], 'Median', 'Mean');
l.Position = [0.05 0.77 0.07 0.09];
l.Box = 'off';

subplot(1,4,2)
[ax2] = fplot_boxplot_log(M.ratio_resil_corr(~isnan(M.growth)), M.growth(~isnan(M.growth)));
ax2.YLabel.String = 'Resilience_{intr.}';
ax2.YLim = [-3.2 6.2];

subplot(1,4,3)
[ax3] = fplot_boxplot_log(M.ratio_resist_corr(~isnan(M.growth)), M.growth(~isnan(M.growth)));
ax3.YLabel.String = 'Resistance_{intr.}';
ax3.YLim = [-0.1 2.5];

subplot(1,4,4)
[ax4] = fplot_boxplot_log(M.ratio_temp_stab_corr(~isnan(M.growth)), M.growth(~isnan(M.growth)));
ax4.YLabel.String = 'Temporal stability_{intr.}';
ax4.YLim = [0.3 50];
ax4.YScale = 'log';

annotation('textbox', [0.045, 0.77, 0, 0], 'string', 'a', 'FontWeight', 'bold')
annotation('textbox', [0.265, 0.77, 0, 0], 'string', 'b', 'FontWeight', 'bold')
annotation('textbox', [0.485, 0.77, 0, 0], 'string', 'c', 'FontWeight', 'bold')
annotation('textbox', [0.705, 0.77, 0, 0], 'string', 'd', 'FontWeight', 'bold')

set(findall(gcf,'-property','FontSize'),'FontSize',16)

n_recov = sum(~isnan(M.ratio_recov_corr(~isnan(M.growth))));
n_resil = sum(~isnan(M.ratio_resil_corr(~isnan(M.growth))));
n_resist = sum(~isnan(M.ratio_resist_corr(~isnan(M.growth))));
n_temp_stab = sum(~isnan(M.ratio_temp_stab_corr(~isnan(M.growth))));

annotation('textbox', [0.155, 0.5, 0.2, 0.2], 'string', ...
    sprintf('n = %d', n_recov), 'Fontsize', 12, 'linestyle', 'none')
annotation('textbox', [0.375, 0.5, 0.2, 0.2], 'string', ...
    sprintf('n = %d', n_resil), 'Fontsize', 12, 'linestyle', 'none')
annotation('textbox', [0.375+0.22, 0.5, 0.2, 0.2], 'string', ...
    sprintf('n = %d', n_resist), 'Fontsize', 12, 'linestyle', 'none')
annotation('textbox', [0.375+0.22*2, 0.5, 0.2, 0.2], 'string', ...
    sprintf('n = %d', n_temp_stab), 'Fontsize', 12, 'linestyle', 'none')

% add organism group labels
xpos = unique(M.growth(~isnan(M.growth)));
for o = 1:length(xpos)
    grpname = unique(M.organism(M.growth==xpos(o)));
    ht = text(ax4, xpos(o), ax4.YLim(2)*1.2, grpname, 'FontSize', 12);
    set(ht,'Rotation',45)
end

ax1.Position(1) = 0.05;
ax2.Position(1) = ax1.Position(1) + 0.22;
ax3.Position(1) = ax1.Position(1) + 0.22*2;
ax4.Position(1) = ax1.Position(1) + 0.22*3;
ax1.Position(2) = 0.25;
ax2.Position(2) = 0.25;
ax3.Position(2) = 0.25;
ax4.Position(2) = 0.25;
ax1.Position(3) = 0.15;
ax2.Position(3) = 0.15;
ax3.Position(3) = 0.15;
ax4.Position(3) = 0.15;
ax1.Position(4) = 0.45;
ax2.Position(4) = 0.45;
ax3.Position(4) = 0.45;
ax4.Position(4) = 0.45;

ax1.XLabel.Position(2) = -30;
ax2.XLabel.Position(2) = -30;
ax3.XLabel.Position(2) = -30;
ax4.XLabel.Position(2) = -30;


%%  Figure S2: Boxplots of published stability
% These are the stability measures as reported in the Hillebrand2020 study.
% In the paper, we use self-calculated estimates based on log-response
% ratios.

figure('color', 'white', 'position', [45,555,1146,380])
subplot(1,4,1)
[ax1] = fplot_boxplot_log(M.recov_LRR(~isnan(M.growth)), M.growth(~isnan(M.growth)));
ax1.YLabel.String = 'Recovery';
ax1.YLim = [-1.5 1.5];
lines =  findobj(gca,'Type','line');
l = legend(ax1, [lines(1) lines(2)], 'Median', 'Mean');
l.Position = [0.05 0.77 0.07 0.09];
l.Box = 'off';

subplot(1,4,2)
[ax2] = fplot_boxplot_log(M.resil_rma_day(~isnan(M.growth)), M.growth(~isnan(M.growth)));
ax2.YLabel.String = 'Resilience';
ax2.YLim = [-1 2];

subplot(1,4,3)
[ax3] = fplot_boxplot_log(M.resist_LRR(~isnan(M.growth)), M.growth(~isnan(M.growth)));
ax3.YLabel.String = 'Resistance';
ax3.YLim = [-4 1];

subplot(1,4,4)
[ax4] = fplot_boxplot_log(M.temp_stab_rma_day(~isnan(M.growth)), M.growth(~isnan(M.growth)));
ax4.YLabel.String = 'Temporal stability';
ax4.YLim = [0.5 300];
ax4.YScale = 'log';

annotation('textbox', [0.045, 0.77, 0, 0], 'string', 'a', 'FontWeight', 'bold')
annotation('textbox', [0.265, 0.77, 0, 0], 'string', 'b', 'FontWeight', 'bold')
annotation('textbox', [0.485, 0.77, 0, 0], 'string', 'c', 'FontWeight', 'bold')
annotation('textbox', [0.705, 0.77, 0, 0], 'string', 'd', 'FontWeight', 'bold')

set(findall(gcf,'-property','FontSize'),'FontSize',16)

% add organism group labels
xpos = unique(M.growth(~isnan(M.growth)));
for o = 1:length(xpos)
    grpname = unique(M.organism(M.growth==xpos(o)));
    ht = text(ax4, xpos(o), ax4.YLim(2)*1.45, grpname, 'FontSize', 12);
    set(ht,'Rotation',45)
end

ax1.Position(1) = 0.05;
ax2.Position(1) = ax1.Position(1) + 0.22;
ax3.Position(1) = ax1.Position(1) + 0.22*2;
ax4.Position(1) = ax1.Position(1) + 0.22*3;
ax1.Position(2) = 0.25;
ax2.Position(2) = 0.25;
ax3.Position(2) = 0.25;
ax4.Position(2) = 0.25;
ax1.Position(3) = 0.15;
ax2.Position(3) = 0.15;
ax3.Position(3) = 0.15;
ax4.Position(3) = 0.15;
ax1.Position(4) = 0.45;
ax2.Position(4) = 0.45;
ax3.Position(4) = 0.45;
ax4.Position(4) = 0.45;

ax1.XLabel.Position(2) = -30;
ax2.XLabel.Position(2) = -30;
ax3.XLabel.Position(2) = -30;
ax4.XLabel.Position(2) = -30;


%% Figure S3: Differences between fast-growers and slow-growers

fasties = ismember(M.organism, {'microbes', 'phytoplankton'});
slowies = ismember(M.organism, {'macrophytes', 'grasses'});

fc = [.5 .5 .5];
ec = 'none';
bw = 0.5;
do = 'ascend';
or = 'horizontal';

figure('color', 'white', 'position', [394,774,1031,204])
subplot(2,5,1)
histogram(categorical(M.habitat(slowies)), 'Facecolor', fc, 'Edgecolor', ec, 'barwidth', bw, 'displayorder', do, 'Orientation', or)
ylabel([{'\bfSlow-growing'};{'group'}])
title('Habitat', 'Fontsize', 10)
ax1 = gca;
subplot(2,5,2)
histogram(categorical(M.dist_group(slowies)), 'Facecolor', fc, 'Edgecolor', ec, 'barwidth', bw, 'displayorder', do, 'Orientation', or)
title('Dist. type', 'Fontsize', 10)
ax2 = gca;
subplot(2,5,3)
histogram(categorical(M.resp(slowies)), 'Facecolor', fc, 'Edgecolor', ec, 'barwidth', bw, 'displayorder', do, 'Orientation', or)
title('Response', 'Fontsize', 10)
ax3 = gca;
subplot(2,5,4)
histogram(categorical(M.open(slowies)), 'Facecolor', fc, 'Edgecolor', ec, 'barwidth', bw, 'displayorder', do, 'Orientation', or)
title('Dispersal', 'Fontsize', 10)
ax4 = gca;

linkaxes([ax1 ax2 ax3 ax4], 'x')

subplot(2,5,6)
histogram(categorical(M.habitat(fasties)), 'Facecolor', fc, 'Edgecolor', ec, 'barwidth', bw, 'displayorder', do, 'Orientation', or)
ylabel([{'\bfFast-growing'};{'group'}])
xlabel('n')
ax6 = gca;
subplot(2,5,7)
histogram(categorical(M.dist_group(fasties)), 'Facecolor', fc, 'Edgecolor', ec, 'barwidth', bw, 'displayorder', do, 'Orientation', or)
xlabel('n')
ax7 = gca;
subplot(2,5,8)
histogram(categorical(M.resp(fasties)), 'Facecolor', fc, 'Edgecolor', ec, 'barwidth', bw, 'displayorder', do, 'Orientation', or)
xlabel('n')
ax8 = gca;
subplot(2,5,9)
histogram(categorical(M.open(fasties)), 'Facecolor', fc, 'Edgecolor', ec, 'barwidth', bw, 'displayorder', do, 'Orientation', or)
xlabel('n')
ax9 = gca;

linkaxes([ax6 ax7 ax8 ax9], 'x')

annotation('textbox', [0.182, 1, 0, 0], 'string', 'a', 'FontWeight', 'bold')
annotation('textbox', [0.345, 1, 0, 0], 'string', 'b', 'FontWeight', 'bold')
annotation('textbox', [0.508, 1, 0, 0], 'string', 'c', 'FontWeight', 'bold')
annotation('textbox', [0.67, 1, 0, 0], 'string', 'd', 'FontWeight', 'bold')
annotation('textbox', [0.182, 0.53, 0, 0], 'string', 'e', 'FontWeight', 'bold')
annotation('textbox', [0.345, 0.53, 0, 0], 'string', 'f', 'FontWeight', 'bold')
annotation('textbox', [0.508, 0.53, 0, 0], 'string', 'g', 'FontWeight', 'bold')
annotation('textbox', [0.67,  0.53, 0, 0], 'string', 'h', 'FontWeight', 'bold')

export_fig(sprintf('%s/differences_slowies_to_fasties.png', fdir), '-r150')

fprintf('\n\nMean latitude for slow growers is %1.0f, and for fast growers %1.0f.',...
    mean(M.lat(slowies)), mean(M.lat(fasties)))


