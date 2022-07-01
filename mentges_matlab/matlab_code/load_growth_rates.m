%% Adds the compiled growth rate estimates
%
% For each organism group included in the Hillebrand & Kunze meta-data (T),
% growth rate estimates are loaded from a compilation of published growth
% rates (compiled by Andrea Mentges 2020).
% References to original publications are given in the table G.
%
% INPUT: 
% - M, containing the meta-data of the Hillebrand and Kunze 2020
%   data set with observations in rows and meta-data in columns
%
% OUTPUT:
% - M, with additional column "growth" (mean of growth rate estimates for
% organism group) and "growth_std" (standard deviation of growth rates for
% the organism group)
% - Table_S3, corresponding to table in supplement, comprising the
% individual growth rate estimates that were used to calculate mean growth
% rate for each organism class

function [M, Table_S3] = load_growth_rates(M, varargin)
%% Load growth rate for organism groups from compilation spreadsheet

% load compilation table
Gorg = readtable('compiled_growth_rates.xlsx');
nrows = size(Gorg,1);

% restrict to estimates with the unit 1/d
% add estimates to new column, either the orginal or, if available,
% the converted estimate (converted to 1/d)
Gorg.estimate = NaN(nrows,1);
Gorg.estimate(strcmp(Gorg.original_unit, '1/d')) = Gorg.original_estimate(strcmp(Gorg.original_unit, '1/d'));
Gorg.estimate(strcmp(Gorg.converted_unit, '1/d')) = Gorg.converted_estimate(strcmp(Gorg.converted_unit, '1/d'));

% remove all estimates whose relation to instantaneous growth rate (IGR)
% is "unclear"
Gorg(strcmp(Gorg.unified_measure, 'unclear'),:) = [];

% restrict table to rows with valid estimates (non-NaN)
G = Gorg;
G(isnan(G.estimate),:) = [];

% get organism order (decreasing growth rate)
[means, ~, ~, names] = grpstats(G.estimate, G.organism);
[B,I] = sort(means, 'descend');

if ~isempty(varargin)
    
    if strcmp(varargin, 'plot')
        % Plot mean growth rates
        figure()
        boxplot(G.estimate, G.organism, 'GroupOrder', names(I))
        set(gca, 'YScale', 'log')
        ylabel('Growth rate [1/d]')
    end
    
end

%% Add the mean growth rates of each group to Hillebrand2020 table

% assign empty columns for mean and std of estimates for each organism
M.growth     = NaN(size(M.organism));
M.growth_std = NaN(size(M.organism));
organisms = unique(G.organism);

for o = 1:length(organisms)
    organism = organisms(o);
    estimates = G.estimate(strcmp(G.organism, organism));
    M.growth(strcmp(M.organism, organism)) = mean(estimates);
    M.growth_std(strcmp(M.organism, organism)) = std(estimates);
end    

% find studies per organism group used for the analyses
used_organisms = unique(M.organism(~isnan(M.growth)));
for u = 1:length(used_organisms)
    organism = used_organisms(u);
    references = unique(G.reference(strcmp(G.organism, organism)));
end   

% set up table with the used growth rates
G_used = G(ismember(G.organism, unique(M.organism)),:);
G_used.Notes = [];
G_used.n_treatment = [];
G_used.scope = [];
G_used.converted_unit = [];
G_used.estimate = [];
G_used.unified_measure = [];
Table_S3 = G_used;

end