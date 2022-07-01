%% Helper function for plotting of boxplots
% Plots the values in y as boxplots for each group along log(x), and
% adds two colored lines indicating mean and median.

function [axh] = fplot_boxplot_log(y,x, varargin)

xscale = 'log';
if ~isempty(varargin)
    xscale = varargin{:};
end

boxplot(y, x,...
    'positions', x, 'widths', 0.01,...
    'PlotStyle', 'compact', 'Colors', [0 0 0], 'MedianStyle', 'line',...
    'OutlierSize', 2, 'Symbol', '.', 'Jitter', 0, 'BoxStyle', 'filled',...
    'Labels', unique(x));
axh = gca;
axh.XLim = [0.063 1.7];
axh.XScale = xscale;
axh.XTick = 10.^[-1 0];
axh.XTickLabels = cellstr(num2str(round(log10(axh.XTick(:))), '10^{%d}'));
axh.XLabel.String = 'Growth rate [1/d]';
hold on
[means, ~, ~, names] = grpstats(y,x);
[medians] = grpstats(y,x, 'median');
x = unique(x);
plot(x(~isnan(means)), means(~isnan(means)), 'color', mycolors('redorange'), 'linewidth', 1.5)
plot(x(~isnan(medians)), medians(~isnan(medians)), 'color', mycolors('yellow'), 'linewidth', 1.5)

end