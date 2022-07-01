%% Calcualtes stability estimates for Hillebrand & Kunze data
%
% Loads the time series of log-response ratios from the orginial study
% (available at datadryad.org/stash/dataset/doi:10.5061/dryad.cz8w9gj09)
% and calculates raw and "corrected" stability estimates from this.
%
% OUTPUT:
% - structure "M" with additional columns for the raw and corrected
% stability estimates (based on the log-response ratios)
% - Table "Table_S2", corresponding to the table in the supplement,
% comprising for each of the seven organism groups from the paper the mean
% and standard deviation of the growth rate estimates
%
% Note: Day 0 is before the disturbance.

function [M, Table_S2] = load_stability_data(M)
%% Load individual time points treatment/control ratios

% import the time series of log-response ratios data
S = readtable('all effect sizes.csv');

% ratio of treatment to control (<1 means treatment has a negative effect,
% i.e. disturbance decreased the biomass)
S.ratio = exp(S.LRR);

%% Derive raw stability measures

% pre-assign columns in the master table for the raw estimates
M.ratio_resil = NaN(size(M.growth));
M.ratio_resil_regression_log = NaN(size(M.growth));
M.ratio_resil_regression_lin = NaN(size(M.growth));
M.ratio_resist = NaN(size(M.growth));
M.ratio_recov = NaN(size(M.growth));
M.ratio_temp_stab = NaN(size(M.growth));

% pre-assign columns in the master table for the corrected estimates
M.ratio_resil_corr = NaN(size(M.growth));
M.ratio_resist_corr = NaN(size(M.growth));
M.ratio_recov_corr = NaN(size(M.growth));
M.ratio_temp_stab_corr = NaN(size(M.growth));

% pre-assign columns in the master table for to find mean sampling interval
% for resistance re-scaling
M.interval = NaN(size(M.growth));

% pre-assign columns in the ratio table for the rescaled day and ratio
S.DAY_resc     = NaN(size(S,1),1);
S.ratio_resc   = NaN(size(S,1),1);
S.DAY_resc_max = NaN(size(S,1),1); % to find the last joint time point

% loop through studies and responses in R (ratio table)
unique_cases = unique(S.caseID);
for c = 1:length(unique_cases)
    caseID = unique_cases(c);
    isCase = strcmp(S.caseID, caseID);
    unique_responses = unique(S.resp_cat(isCase));
    
    for r = 1:length(unique_responses)
        response = unique_responses(r);
        isResponse = strcmp(S.resp_cat, response);
        day = S.DAY(isCase & isResponse);
        ratio = S.ratio(isCase & isResponse);
        
        % find respective row in T table
        idx = strcmp(M.caseID, caseID) & strcmp(M.resp_cat, response);
        assert(sum(idx)<=1, 'Multiple rows match!')
        
        %% Calculate raw stability
        
        % post-disturbance time points (day 0 is before the disturbance)
        post = day>0;
        
        % get day-to-day slopes
        xdnndt = abs(ratio(post)-1); % center equilibrium around zero
        midxdnndt = (xdnndt(1:(end-1))+xdnndt(2:end))/2;
        slopes = diff(xdnndt)./diff(day(post))./midxdnndt;
        
        % get overall slope (like hillebrand2020)
        regression = []; 
        if sum(post)>1 % if more than 1 time point after disturbance
%            slopes_first_last = (ratio(end)-ratio(1))/(day(end)-day(1));
            x = log(day(post));
            y = ratio(post);
            p = polyfit(x,y,1);
            regression_log = p(1);
            
            x = day(post);
            y = ratio(post);
            p = polyfit(x,y,1);
            regression_lin = p(1);
        end
        
        % resilience = max of day to day linear slopes
        if sum(post)>1 % if more than 1 time point after disturbance
            M.ratio_resil(idx) = max(slopes);
%             T.ratio_resil_first_last(idx) = slopes_first_last;
%             T.ratio_resil_first_last(idx) = mean(slopes);
            M.ratio_resil_regression_log(idx) = regression_log;
            M.ratio_resil_regression_lin(idx) = regression_lin;
        end
        
        % resistance = ratio at first time point after disturbance
        M.ratio_resist(idx) = ratio(find(post==1, 1, 'first'));
       
        % recovery = ratio at last day
        M.ratio_recov(idx) = ratio(find(post==1, 1, 'last'));
        
        % temporal stability = 1/CV of ratio
        if length(day)>1
            M.ratio_temp_stab(idx) = 1/(std(ratio)/mean(ratio));
        end
        
        
        %% Derive re-scaled time axis for resilience and resistance
        
        % get growth rate
        growth_rate = M.growth(idx);

        % if a growth rate exists and at least two time points after
        % disturbance, rescale
        t_dist = min(day(post));
        if ~isempty(growth_rate) && sum(post)>1
            day_resc = (day(post)-t_dist)*growth_rate;
            S.DAY_resc(isCase & isResponse) = [NaN(sum(~post,1));day_resc];
            S.DAY_resc_max(isCase & isResponse) = max(day_resc);
        end
        
        % finding first time point after disturbance for resistance
        % correction
        M.interval(idx) = min(day(post));     
        
    end
end

%% Calculate corrected stability (accounting for growth rate)

% This is done in an extra loop, because first the time points need to be
% re-scaled and the maximum joint time point of all time series needs to be
% identified
last_joint = min(S.DAY_resc_max);

% loop through studies and responses in R (ratio table)
unique_cases = unique(S.caseID);
for c = 1:length(unique_cases)
    
    caseID = unique_cases(c);
    isCase = strcmp(S.caseID, caseID);
    unique_responses = unique(S.resp_cat(isCase));
    
    for r = 1:length(unique_responses)
        
        response = unique_responses(r);
        isResponse = strcmp(S.resp_cat, response);
        day = S.DAY(isCase & isResponse);
        day_resc = S.DAY_resc(isCase & isResponse);
        
        % find respective row in T table
        idx = strcmp(M.caseID, caseID) & strcmp(M.resp_cat, response);
        assert(sum(idx)<=1, 'Multiple rows match!')

        % restrict rescaled time to joint time points
        idx_joint  = (day_resc<=last_joint) & day_resc >= 0;
        day_resc   = day_resc(idx_joint);
        ratio      = S.ratio(isCase & isResponse); 
        ratio_resc = ratio(idx_joint);
        
        % get corrected resilience and recovery
        xdnndt = abs(ratio_resc-1);
        midxdnndt = (xdnndt(1:(end-1))+xdnndt(2:end))/2;
        slopes = diff(xdnndt)./diff(day_resc)./midxdnndt;
        if length(slopes)>=1
            M.ratio_resil_corr(idx) = max(slopes);
        end
        if ~isempty(day_resc)
            M.ratio_recov_corr(idx) = ratio_resc(end);
        end
        
        % get corrected resistance (calculated with original, not rescaled
        % time, but adapted sampling interval)
        growth_rate = M.growth(idx);
        if length(day)>1 & ~isnan(growth_rate)
            % sampling interval is adapted relative to the fastest growing
            % organism
            r_fast = max(M.growth);
            cases_fast = M.growth == r_fast;
            % minimum time passed between first sampling of fastest organism and disturbance
            timestep_fast = min(M.interval(cases_fast));
            % this is the time that should pass between samplings for the slow grower
            timestep_slow = timestep_fast*r_fast/growth_rate;
            % interpolate respective time step from slow grower data
            x = day;
            y = ratio;
            x2 = timestep_slow;
            y2 = interp1(x, y, x2, 'linear');
            M.ratio_resist_corr(idx) = y2;

        end
        
        % get corrected temporal stability
        if length(day)>1
            M.ratio_temp_stab_corr(idx) = 1/((std(ratio)/mean(ratio))*sqrt(growth_rate));
        end

    end
end


%% Number of data points for the different stability metrics

% number of data points in published estimates
sum(~isnan(M.recov_LRR(~isnan(M.growth))));
sum(~isnan(M.resil_rma_day(~isnan(M.growth))));
sum(~isnan(M.resist_LRR(~isnan(M.growth))));
sum(~isnan(M.temp_stab_rma_day(~isnan(M.growth))));


% number of data points in self-calculated raw estimates
sum(~isnan(M.ratio_recov(~isnan(M.growth))));
sum(~isnan(M.ratio_resil(~isnan(M.growth))));
sum(~isnan(M.ratio_resist(~isnan(M.growth))));
sum(~isnan(M.ratio_temp_stab(~isnan(M.growth))));

% number of data points in corrected estimates
sum(~isnan(M.ratio_recov_corr(~isnan(M.growth))));
sum(~isnan(M.ratio_resil_corr(~isnan(M.growth))));
sum(~isnan(M.ratio_resist_corr(~isnan(M.growth))));
sum(~isnan(M.ratio_temp_stab_corr(~isnan(M.growth))));

% Table S1: show used growth rate means
[a,b] = unique(M.growth);
M.organism(b);
[c,d] = unique(M.growth_std);
tableS1 = cell2table(M.organism(b), 'VariableNames', {'Organism_group'});
tableS1.Mean_growth_rate = a;
tableS1.Std_growth_rate = c;
tableS1 = sortrows(tableS1,2);
Table_S2 = tableS1;

%% Data cleaning of ratio estimates: removing outliers

% removing outliers based on resistance
while any(M.ratio_resist/nanmedian(M.ratio_resist)>50)
    fprintf('\nRemoving the following outlier from raw resistance metric, which is %1.2f times higher than median:',...
        nanmax(M.ratio_resist/nanmedian(M.ratio_resist)))
    [~,ind] = nanmax(M.ratio_resist);
    % % alternatively: just put which study to remove
    % ind = find(strcmp(T.caseID, 'HH024_3'));
    M.ref(ind); % this is the respective reference
    M(strcmp(M.ref, M.ref(ind)),:); % these are all the data from that study
    M(ind,:) % this is the row that will be removed
    M(ind,:) = [];
end

% Note: outlier studies are completely excluded from the table, 
% i.e. for all stability measures 

%% Print number of marine, freshwater, and terrestrial samples

n_mari = sum(strcmp(M.system, 'marine'));
n_fres = sum(strcmp(M.system, 'freshwater'));
n_terr = sum(strcmp(M.system, 'terrestrial'));

fprintf('\nOut of the %d total observations, %d were marine, %d were freshwater, and %d were terrestrial',...
    size(M,1), n_mari, n_fres, n_terr)

close all

end